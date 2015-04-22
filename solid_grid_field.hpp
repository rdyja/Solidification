#include <Grid/gridfield.h>
#include "solid_node_data.hpp"

class SolidInputData;

class SolidGridField:public TALYFEMLIB::GridField<SolidNodeData>
{
public:
    SolidGridField(const SolidInputData& inputData) :  inputData_(inputData) {
    }

    void SetIC(int);
    void SetBndrIndicator(double width, double length, double height);
private:
    const SolidInputData& inputData_;
};
