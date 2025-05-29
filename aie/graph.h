// By: Austin Owens
// Date: 6/3/2024
// Desc: Using Vitis DSP lib to perform 1D FFT operation

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
        kernel px_demux_km[1];

        // Image reconstruction kernel module
        kernel img_rec_km[IMG_SOLVERS];

        //***** PACKET SWITCHING OBJECTS *****//
        pktsplit<IMG_SOLVERS> sp;
        pktmerge<IMG_SOLVERS> mg;

    public:
        //***** GMIO PORT OBJECTS *****//

        // Data broadcaster GMIO ports
        input_gmio gmio_in_st;
        input_gmio gmio_in_rc;

        // Pixel demux GMIO ports
        input_gmio gmio_in_xyz_px[1];


        //***** PLIO PORT OBJECTS *****//

        // Packet router PLIO ports
        output_plio plio_pkt_rtr_out[1];


        //***** RTP PORT OBJECTS *****//
        input_port rtp_dump_img_in[IMG_SOLVERS];

        BackProjectionGraph() {

            //***** KERNELS *****//
 
            // Data broadcaster kernel
            data_bc_km = kernel::create(data_broadcast_kern);

            // Pixel demux kernel
            px_demux_km[0] = kernel::create(px_demux_kern);
            
            // Image reconstruct kernel
            for (int i=0; i<IMG_SOLVERS; i++) {
                img_rec_km[i] = kernel::create_object<ImgReconstruct>(i);
            }

            //***** PACKET SWITCHING OBJECTS *****//
            
            // Packet spliter/merger
            sp = pktsplit<IMG_SOLVERS>::create();
            mg = pktmerge<IMG_SOLVERS>::create();


            //***** GMIO PORTS *****//

            // Slow time splicer GMIO ports
            gmio_in_st = input_gmio::create("gmio_in_st_" + std::to_string(bp_graph_insts), 256, 1000);
            gmio_in_rc = input_gmio::create("gmio_in_rc_" + std::to_string(bp_graph_insts), 256, 1000);

            // Pixel demux GMIO ports
            gmio_in_xyz_px[0] = input_gmio::create("gmio_in_xyz_px_" + std::to_string(bp_graph_insts), 256, 1000);
            

            //***** PLIO PORTS *****//
            std::string data_file_str = "aie_to_plio_switch_" + std::to_string(0) + ".csv";
            std::string port_name_str = "plio_pkt_rtr_out_" + std::to_string(bp_graph_insts) + "_" + std::to_string(0);
            plio_pkt_rtr_out[0] = output_plio::create(port_name_str.c_str(), plio_128_bits, data_file_str.c_str());


            //***** GMIO CONNECTIONS *****//

            // GMIO x, y, z and ref range to data broadcaster kernel
            connect(gmio_in_st.out[0], data_bc_km.in[0]);
            connect(gmio_in_rc.out[0], data_bc_km.in[1]);

            // Pixel GMIO ports pixel demux kernel
            connect(gmio_in_xyz_px[0].out[0], px_demux_km[0].in[0]);

            // Packet merger to PLIO packet router
            connect(mg.out[0], plio_pkt_rtr_out[0].in[0]);


            //***** AIE TO AIE CONNECTIONS *****//

            // Slow time splicer to image reconstruction
            for (int i=0; i<IMG_SOLVERS; i++) {
                connect(data_bc_km.out[0], img_rec_km[i].in[0]);
                connect(data_bc_km.out[1], img_rec_km[i].in[1]);
            }

            // Packet splitter to image reconstruction
            for (int i=0; i<IMG_SOLVERS; i++) {
                connect(sp.out[i], img_rec_km[i].in[2]);
            }

            // Pixel demux to packet splitter
            connect(px_demux_km[0].out[0], sp.in[0]);

            // Image reconstruction to packet merger
            for (int i=0; i<IMG_SOLVERS; i++) {
                connect(img_rec_km[i].out[0], mg.in[i]);
            }


            //***** RTP CONNECTIONS *****//
            
            // Image reconstruction to valid_bounds RTP param
            for (int i=0; i<IMG_SOLVERS; i++) {
                connect<parameter>(rtp_dump_img_in[i], img_rec_km[i].in[3]);
            }


            //***** SOURCE FILES *****//

            source(data_bc_km) = "backprojection.cc";
            source(px_demux_km[0]) = "backprojection.cc";
            for (int i=0; i<IMG_SOLVERS; i++)
                source(img_rec_km[i]) = "backprojection.cc";


            //***** RUNTIME RATIOS *****//

            runtime<ratio>(data_bc_km) = 1.0;
            runtime<ratio>(px_demux_km[0]) = 1.0;
            for (int i=0; i<IMG_SOLVERS; i++)
                runtime<ratio>(img_rec_km[i]) = 1.0;

            ++bp_graph_insts;
        
        }
};
