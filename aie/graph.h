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

        // Slow time splicer kernel module
        kernel sts_km;

        // Pixel arbiter kernel module
        kernel px_arb_km[1];

        // Image reconstruction kernel module
        kernel img_rec_km[IMG_SOLVERS];

        //***** PACKET SWITCHING OBJECTS *****//
        pktsplit<IMG_SOLVERS> sp;
        pktmerge<IMG_SOLVERS> mg;

    public:
        //***** GMIO PORT OBJECTS *****//

        // Slow time splicer GMIO ports
        input_gmio gmio_in_x_ant_pos;
        input_gmio gmio_in_y_ant_pos;
        input_gmio gmio_in_z_ant_pos;
        input_gmio gmio_in_ref_range;
        input_gmio gmio_in_rc;

        // Pixel arbiter GMIO ports
        input_gmio gmio_in_xyz_px[1];

        // Image reconstruction GMIO ports
        input_gmio gmio_in_xy_px[IMG_SOLVERS];
        input_gmio gmio_in_z_px[IMG_SOLVERS];
        //output_gmio gmio_out_img[1];


        //***** PLIO PORT OBJECTS *****//

        // Packet router PLIO ports
        output_plio plio_pkt_rtr_out[1];


        //***** RTP PORT OBJECTS *****//
        input_port rtp_dump_img_in[IMG_SOLVERS];

        BackProjectionGraph() {

            //***** KERNELS *****//
 
            // Slow time splicer kernel
            sts_km = kernel::create(slowtime_splicer_kern);

            // Pixel arbiter kernel
            px_arb_km[0] = kernel::create(px_arbiter_kern);
            
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
            gmio_in_x_ant_pos = input_gmio::create("gmio_in_x_ant_pos_" + std::to_string(bp_graph_insts), 256, 1000);
            gmio_in_y_ant_pos = input_gmio::create("gmio_in_y_ant_pos_" + std::to_string(bp_graph_insts), 256, 1000);
            gmio_in_z_ant_pos = input_gmio::create("gmio_in_z_ant_pos_" + std::to_string(bp_graph_insts), 256, 1000);
            gmio_in_ref_range = input_gmio::create("gmio_in_ref_range_" + std::to_string(bp_graph_insts), 256, 1000);
            gmio_in_rc = input_gmio::create("gmio_in_rc_" + std::to_string(bp_graph_insts), 256, 1000);

            // Pixel arbiter GMIO ports
            gmio_in_xyz_px[0] = input_gmio::create("gmio_in_xyz_px_" + std::to_string(bp_graph_insts), 256, 1000);
            
            // Image reconstruct GMIO ports
            //for (int i=0; i<IMG_SOLVERS; i++) {
            //    gmio_out_img[i] = output_gmio::create("gmio_out_img_" + std::to_string(bp_graph_insts) + "_" + std::to_string(i), 256, 1000);
            //}

            // Output image GMIO ports
            //gmio_out_img[0] = output_gmio::create("gmio_out_img_" + std::to_string(bp_graph_insts) + "_" + std::to_string(0), 256, 1000);


            //***** PLIO PORTS *****//
            std::string data_file_str = "aie_to_plio_switch_" + std::to_string(0) + ".csv";
            std::string port_name_str = "plio_pkt_rtr_out_" + std::to_string(bp_graph_insts) + "_" + std::to_string(0);
            plio_pkt_rtr_out[0] = output_plio::create(port_name_str.c_str(), plio_128_bits, data_file_str.c_str());


            //***** GMIO CONNECTIONS *****//

            // GMIO x, y, z and ref range to slow time splicer kernel
            connect(gmio_in_x_ant_pos.out[0], sts_km.in[0]);
            connect(gmio_in_y_ant_pos.out[0], sts_km.in[1]);
            connect(gmio_in_z_ant_pos.out[0], sts_km.in[2]);
            connect(gmio_in_ref_range.out[0], sts_km.in[3]);
            connect(gmio_in_rc.out[0], sts_km.in[4]);

            // Pixel GMIO ports pixel arbiter kernel
            connect(gmio_in_xyz_px[0].out[0], px_arb_km[0].in[0]);

            // Packet merger to GMIO output image
            //connect(mg.out[0], gmio_out_img[0].in[0]);

            // Packet merger to PLIO packet router
            connect(mg.out[0], plio_pkt_rtr_out[0].in[0]);

            // Image reconstruction kernel output to GMIO
            //for (int i=0; i<IMG_SOLVERS; i++) {
            //    connect(img_rec_km[i].out[0], gmio_out_img[i].in[0]);
            //}


            //***** AIE TO AIE CONNECTIONS *****//

            // Slow time splicer to image reconstruction
            for (int i=0; i<IMG_SOLVERS; i++) {
                connect(sts_km.out[0], img_rec_km[i].in[0]);
                connect(sts_km.out[1], img_rec_km[i].in[1]);
            }

            // Packet splitter to image reconstruction
            for (int i=0; i<IMG_SOLVERS; i++) {
                connect(sp.out[i], img_rec_km[i].in[2]);
            }

            // Pixel arbiter to packet splitter
            connect(px_arb_km[0].out[0], sp.in[0]);

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

            source(sts_km) = "backprojection.cc";
            source(px_arb_km[0]) = "backprojection.cc";
            for (int i=0; i<IMG_SOLVERS; i++)
                source(img_rec_km[i]) = "backprojection.cc";


            //***** RUNTIME RATIOS *****//

            runtime<ratio>(sts_km) = 1.0;
            runtime<ratio>(px_arb_km[0]) = 1.0;
            for (int i=0; i<IMG_SOLVERS; i++)
                runtime<ratio>(img_rec_km[i]) = 1.0;

            ++bp_graph_insts;
        
        }
};
