// By: Austin Owens
// Date: 6/3/2024
// Desc: ADF graph for backprojection

#pragma once

#include "adf.h"
#include "custom_kernels.h"

using namespace adf;

extern uint8_t bp_graph_insts;
uint8_t bp_subgraph_insts;

class BackProjectionSubgraph: public graph {
    public:

        //***** KERNEL OBJECTS *****//

        // Pixel demux kernel module
        kernel px_demux_km;

        // Image reconstruction kernel module
        kernel img_rec_km[IMG_SOLVERS_PER_SWITCH];


        //***** PACKET SWITCHING OBJECTS *****//
        pktsplit<IMG_SOLVERS_PER_SWITCH> sp;
        pktmerge<IMG_SOLVERS_PER_SWITCH> mg;


        //***** GMIO PORT OBJECTS *****//

        // Pixel demux GMIO port
        input_gmio gmio_in_xyz_px;


        //***** PLIO PORT OBJECTS *****//

        // Packet router PLIO port
        output_plio plio_pkt_rtr_out;


        BackProjectionSubgraph() {

            //***** KERNELS *****//
 
            // Pixel demux kernel
            px_demux_km = kernel::create(px_demux_kern);
            
            // Image reconstruct kernels
            for (int i=0; i<IMG_SOLVERS_PER_SWITCH; i++)
                img_rec_km[i] = kernel::create_object<ImgReconstruct>(IMG_SOLVERS_PER_SWITCH*bp_subgraph_insts + i);

            //***** PACKET SWITCHING OBJECTS *****//
            
            // Packet spliter/merger
            sp = pktsplit<IMG_SOLVERS_PER_SWITCH>::create();
            mg = pktmerge<IMG_SOLVERS_PER_SWITCH>::create();


            //***** GMIO PORTS *****//

            // Pixel demux GMIO ports
            std::string xyz_px_str = "gmio_in_xyz_px_" + std::to_string(bp_graph_insts) + "_" + std::to_string(bp_subgraph_insts);
            gmio_in_xyz_px = input_gmio::create(xyz_px_str.c_str(), 256, 1000);
            

            //***** PLIO PORTS *****//

            std::string plio_data_file_str = "aie_to_plio_switch_" + std::to_string(bp_graph_insts) + "_" + std::to_string(bp_subgraph_insts) + ".csv";
            std::string plio_pkt_rtr_str = "plio_pkt_rtr_out_" + std::to_string(bp_graph_insts) + "_" + std::to_string(bp_subgraph_insts);
            plio_pkt_rtr_out = output_plio::create(plio_pkt_rtr_str.c_str(), plio_128_bits, plio_data_file_str.c_str());


            //***** GMIO CONNECTIONS *****//

            // Pixel GMIO ports pixel demux kernel
            connect(gmio_in_xyz_px.out[0], px_demux_km.in[0]);

            // Packet merger to PLIO packet router
            connect(mg.out[0], plio_pkt_rtr_out.in[0]);


            //***** AIE TO AIE CONNECTIONS *****//

            // Packet splitter to image reconstruction
            for (int i=0; i<IMG_SOLVERS_PER_SWITCH; i++)
                connect(sp.out[i], img_rec_km[i].in[2]);

            // Pixel demux to packet splitter
            connect(px_demux_km.out[0], sp.in[0]);

            // Image reconstruction to packet merger
            for (int i=0; i<IMG_SOLVERS_PER_SWITCH; i++) {
                connect(img_rec_km[i].out[0], mg.in[i]);
            }

            // ***** LOCATION CONSTRAINTS ***** //

            // On VC1902 only columns 6-44 are PLIO-capable
            int base_col = 0;
            int col_start = base_col + bp_subgraph_insts * 7;
            int col_end   = col_start + 6;
            int row_start = 0;
            int row_end   = 7;
     
            location<graph>(*this) = area_group({
                { aie_tile, col_start, row_start, col_end, row_end },
                { shim_tile, col_start, 0, col_end, 0 } // Needed for PLIOs
            });

            //***** SOURCE FILES *****//

            source(px_demux_km) = "backprojection.cc";
            for (int i=0; i<IMG_SOLVERS_PER_SWITCH; i++)
                source(img_rec_km[i]) = "backprojection.cc";


            //***** RUNTIME RATIOS *****//

            runtime<ratio>(px_demux_km) = 1.0;
            for (int i=0; i<IMG_SOLVERS_PER_SWITCH; i++)
                runtime<ratio>(img_rec_km[i]) = 1.0;

            bp_subgraph_insts++;
        
        }
};


class BackProjectionGraph: public graph {
    private:
        //***** KERNEL OBJECTS *****//

        // Data broadcaster kernel module
        kernel data_bc_km;
        

    public:
        //***** GRAPH OBJECTS *****//

        // Create multiple subgraphs of backprojection clusters
        BackProjectionSubgraph bpCluster[AIE_SWITCHES];

        //***** GMIO PORT OBJECTS *****//

        // Data broadcaster GMIO ports
        input_gmio gmio_in_st;
        input_gmio gmio_in_rc;

        //***** RTP PORT OBJECTS *****//
        input_port rtp_dump_img_in[IMG_SOLVERS];

        BackProjectionGraph() {

            //***** KERNELS *****//
 
            // Data broadcaster kernel
            data_bc_km = kernel::create(data_broadcast_kern);

            //***** GMIO PORTS *****//

            // Slow time splicer GMIO ports
            gmio_in_st = input_gmio::create("gmio_in_st_" + std::to_string(bp_graph_insts), 256, 1000);
            gmio_in_rc = input_gmio::create("gmio_in_rc_" + std::to_string(bp_graph_insts), 256, 1000);

            //***** GMIO CONNECTIONS *****//

            // GMIO x, y, z and ref range to data broadcaster kernel
            connect(gmio_in_st.out[0], data_bc_km.in[0]);
            connect(gmio_in_rc.out[0], data_bc_km.in[1]);

            //***** AIE TO AIE CONNECTIONS *****//

            // Data broadcaster to image reconstruction
            for (int j=0; j<AIE_SWITCHES; j++) {
                for (int i=0; i<IMG_SOLVERS_PER_SWITCH; i++) {
                    connect(data_bc_km.out[0], bpCluster[j].img_rec_km[i].in[0]);
                    connect(data_bc_km.out[1], bpCluster[j].img_rec_km[i].in[1]);
                    single_buffer(bpCluster[j].img_rec_km[i].in[0]);
                    single_buffer(bpCluster[j].img_rec_km[i].in[1]);
                }
            }

            //***** RTP CONNECTIONS *****//
            
            // Image reconstruction to valid_bounds RTP param
            for (int j=0; j<AIE_SWITCHES; j++) {
                for (int i=0; i<IMG_SOLVERS_PER_SWITCH; i++) {
                    connect<parameter>(rtp_dump_img_in[IMG_SOLVERS_PER_SWITCH*j + i], bpCluster[j].img_rec_km[i].in[3]);
                }
            }

            //***** SOURCE FILES *****//

            source(data_bc_km) = "backprojection.cc";

            //***** RUNTIME RATIOS *****//

            runtime<ratio>(data_bc_km) = 1.0;

            bp_graph_insts++;
        
        }
};

