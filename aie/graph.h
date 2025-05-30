// By: Austin Owens
// Date: 6/3/2024
// Desc: ADF graph for backprojection

#pragma once

#include "adf.h"
#include "custom_kernels.h"

using namespace adf;

extern uint8_t bp_graph_insts;

class BackProjectionGraph: public graph
{
    private:

        //***** KERNEL OBJECTS *****//

        // Data broadcaster kernel module
        kernel data_bc_km;

        // Pixel demux kernel module
        kernel px_demux_km[AIE_SWITCHES];

        // Image reconstruction kernel module
        kernel img_rec_km[IMG_SOLVERS];

        //***** PACKET SWITCHING OBJECTS *****//
        pktsplit<IMG_SOLVERS_PER_SWITCH> sp[AIE_SWITCHES];
        pktmerge<IMG_SOLVERS_PER_SWITCH> mg[AIE_SWITCHES];

    public:
        //***** GMIO PORT OBJECTS *****//

        // Data broadcaster GMIO ports
        input_gmio gmio_in_st;
        input_gmio gmio_in_rc;

        // Pixel demux GMIO ports
        input_gmio gmio_in_xyz_px[AIE_SWITCHES];


        //***** PLIO PORT OBJECTS *****//

        // Packet router PLIO ports
        output_plio plio_pkt_rtr_out[AIE_SWITCHES];


        //***** RTP PORT OBJECTS *****//
        input_port rtp_dump_img_in[IMG_SOLVERS];

        BackProjectionGraph() {

            //***** KERNELS *****//
 
            // Data broadcaster kernel
            data_bc_km = kernel::create(data_broadcast_kern);

            // Pixel demux kernels
            for (int i=0; i<AIE_SWITCHES; i++) {
                px_demux_km[i] = kernel::create(px_demux_kern);
            }
            
            // Image reconstruct kernel
            for (int i=0; i<IMG_SOLVERS; i++) {
                img_rec_km[i] = kernel::create_object<ImgReconstruct>(i);
            }

            //***** PACKET SWITCHING OBJECTS *****//
            
            // Packet spliters/mergers
            for (int i=0; i<AIE_SWITCHES; i++) {
                sp[i] = pktsplit<IMG_SOLVERS_PER_SWITCH>::create();
                mg[i] = pktmerge<IMG_SOLVERS_PER_SWITCH>::create();
            }


            //***** GMIO PORTS *****//

            // Slow time splicer GMIO ports
            gmio_in_st = input_gmio::create("gmio_in_st_" + std::to_string(bp_graph_insts), 256, 1000);
            gmio_in_rc = input_gmio::create("gmio_in_rc_" + std::to_string(bp_graph_insts), 256, 1000);

            // Pixel demux GMIO ports
            for (int i=0; i<AIE_SWITCHES; i++) {
                std::string xyz_px_str = "gmio_in_xyz_px_" + std::to_string(bp_graph_insts) + "_" + std::to_string(i);
                gmio_in_xyz_px[i] = input_gmio::create(xyz_px_str.c_str(), 256, 1000);
            }
            

            //***** PLIO PORTS *****//
            for (int i=0; i<AIE_SWITCHES; i++) {
                std::string data_file_str = "aie_to_plio_switch_" + std::to_string(i) + ".csv";
                std::string pkt_rtr_str = "plio_pkt_rtr_out_" + std::to_string(bp_graph_insts) + "_" + std::to_string(i);
                plio_pkt_rtr_out[i] = output_plio::create(pkt_rtr_str.c_str(), plio_128_bits, data_file_str.c_str());
            }


            //***** GMIO CONNECTIONS *****//

            // GMIO x, y, z and ref range to data broadcaster kernel
            connect(gmio_in_st.out[0], data_bc_km.in[0]);
            connect(gmio_in_rc.out[0], data_bc_km.in[1]);

            // Pixel GMIO ports pixel demux kernel
            for (int i=0; i<AIE_SWITCHES; i++) {
                connect(gmio_in_xyz_px[i].out[0], px_demux_km[i].in[0]);
            }

            // Packet merger to PLIO packet router
            for (int i=0; i<AIE_SWITCHES; i++) {
                connect(mg[i].out[0], plio_pkt_rtr_out[i].in[0]);
            }


            //***** AIE TO AIE CONNECTIONS *****//

            // Slow time splicer to image reconstruction
            for (int i=0; i<IMG_SOLVERS; i++) {
                connect(data_bc_km.out[0], img_rec_km[i].in[0]);
                connect(data_bc_km.out[1], img_rec_km[i].in[1]);
            }

            // Packet splitters to image reconstruction
            for (int j=0; j<AIE_SWITCHES; j++) {
                for (int i=0; i<IMG_SOLVERS_PER_SWITCH; i++) {
                    connect(sp[j].out[i], img_rec_km[IMG_SOLVERS_PER_SWITCH*j + i].in[2]);
                }
            }

            // Pixel demux to packet splitters
            for (int i=0; i<AIE_SWITCHES; i++) {
                connect(px_demux_km[i].out[0], sp[i].in[0]);
            }

            // Image reconstruction to packet mergers
            for (int j=0; j<AIE_SWITCHES; j++) {
                for (int i=0; i<IMG_SOLVERS_PER_SWITCH; i++) {
                    connect(img_rec_km[IMG_SOLVERS_PER_SWITCH*j + i].out[0], mg[j].in[i]);
                }
            }


            //***** RTP CONNECTIONS *****//
            
            // Image reconstruction to valid_bounds RTP param
            for (int i=0; i<IMG_SOLVERS; i++) {
                connect<parameter>(rtp_dump_img_in[i], img_rec_km[i].in[3]);
            }


            //***** SOURCE FILES *****//

            source(data_bc_km) = "backprojection.cc";
            for (int i=0; i<AIE_SWITCHES; i++) {
                source(px_demux_km[i]) = "backprojection.cc";
            }
            for (int i=0; i<IMG_SOLVERS; i++)
                source(img_rec_km[i]) = "backprojection.cc";


            //***** RUNTIME RATIOS *****//

            runtime<ratio>(data_bc_km) = 1.0;
            for (int i=0; i<AIE_SWITCHES; i++) {
                runtime<ratio>(px_demux_km[i]) = 1.0;
            }
            for (int i=0; i<IMG_SOLVERS; i++)
                runtime<ratio>(img_rec_km[i]) = 1.0;

            ++bp_graph_insts;
        
        }
};
