#include <Grid/gridfield.h>
#include "solid_node_data.hpp"
#include "contact_bounds.hpp"

class SolidInputData;

class SolidGridField:public TALYFEMLIB::GridField<SolidNodeData>
{
public:
    SolidGridField(const SolidInputData& inputData) :  inputData_(inputData) {
    }

    void SetIC(int);

private:
    const SolidInputData& inputData_;
};
