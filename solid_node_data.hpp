#ifndef SOLID_NODE_DATA_HPP
#define SOLID_NODE_DATA_HPP

#include <Grid/nodedata.h>

class SolidNodeData:public TALYFEMLIB::NODEData {
	public:
            virtual const double& value(int index) const {
                switch(index) {
                    case 0:
                        return t0;
                    /*case 1:
                        return t1;
                    case 2:
                        return t2;*/
                    case 1:
                        return v;
                    case 2:
                        return fs;
                    case 3:
                        return ts;
                    case 4:
                        return rz;
                    /*case 7:
                        return flux_x;
                    case 8:
                        return flux_y;
                    case 9:
                        return flux_z;*/
                    case 5:
                        return twe;
                    case 6:
                        return capprox;
                }
                throw std::string("Wrong node value");
            }
            virtual double& value(int index) {
                switch(index) {
                    case 0:
                        return t0;
                    /*case 1:
                        return t1;
                    case 2:
                        return t2;*/
                    case 1:
                        return v;
                    case 2:
                        return fs;
                    case 3:
                        return ts;
                    case 4:
                        return rz;
                    /*case 7:
                        return flux_x;
                    case 8:
                        return flux_y;
                    case 9:
                        return flux_z;*/
                    case 5:
                        return twe;
                    case 6:
                        return capprox;
                }
                throw std::string("Wrong node value");
            }
            static const char* name(int index) {
                switch(index) {
                    case 0:
                        return "T";
                    /*case 1:
                        return "t1";
                    case 2:
                        return "t-1";*/
                    case 1:
                        return "v";
                    case 2:
                        return "fs";
                    case 3:
                        return "ts";
                    case 4:
                        return "rz";
                    /*case 7:
                        return "fluxX";
                    case 8:
                        return "fluxY";
                    case 9:
                        return "fluxZ";*/
                    case 5:
                        return "T_with_eutectic";
                    case 6:
                        return "c_approx";
                }
                throw std::string("Wrong node value");
            }
	   double get_prev_temp() const {
               return t0;
           }

           double get_prev_minus_1_temp() const {
               return t2;
           }

           void set_prev_minus_1_temp(double newVal) {
               t2 = newVal;
           }

           void set_curr_temp(double newVal) {
               t0 = newVal;
           }

           void set_velocity(double newVal) {
               v = newVal;
           }
           double get_velocity() const {
               return v;
           }

           double get_solid_fraction() const {
               return fs;
           }

           void set_solid_fraction(double newVal) {
               fs = newVal;
           }

           double get_real_solidus_temperature() const {
               return ts;
           }

           void set_real_solidus_temperature(double newVal) {
               ts = newVal;
           }

           double get_grain_size() const {
               return rz;
           }

           void set_grain_size(double newVal) {
               rz = newVal;
           }

           double get_t_with_eutectic() const {
               return twe;
           }

           void set_t_with_eutectic(double newVal) {
               twe = newVal;
           }

           double get_capprox() const {
               return capprox;
           }

           void set_capprox(double newVal) {
               capprox = newVal;
           }

           static int valueno() {
               return 7;
           }
           void UpdateDataStructures() {
               t2 = t1;
               t1 = t0;
           }

        private:
            double t0, t1, t2, twe, v, fs, ts, rz, flux_x, flux_y, flux_z, capprox;
};

#endif
