#ifndef SOLID_NODE_DATA_HPP
#define SOLID_NODE_DATA_HPP

#include <Grid/nodedata.h>

class SolidNodeData:public TALYFEMLIB::NODEData {
	public:
            virtual double& value(int index) {
                switch(index) {
                    case 0:
                        return t0;
                    case 1:
                        return t1;
                    case 2:
                        return t2;
                    case 3:
                        return v;
                }
                throw std::string("Wrong node value");
            }
            static const char* name(int index) {
                switch(index) {
                    case 0:
                        return "t0";
                    case 1:
                        return "t1";
                    case 2:
                        return "t-1";
                    case 3:
                        return "v";
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
           static int valueno() {
               return 4;
           }
           void UpdateDataStructures() {
               t2 = t1;
               t1 = t0;
           }

        private:
            double t0, t1, t2, v;
};

#endif
