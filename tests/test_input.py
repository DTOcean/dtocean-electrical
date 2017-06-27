
from dtocean_electrical.inputs import ElectricalArrayData


def test_ElectricalArrayData_set_histogram_edges_scalar(mocker):
    
    # Fake a ElectricalMachineData object
    machine_data = mocker.Mock()
    machine_data.power = 1e6
    machine_data.power_factor = 1.
    
    mocker.patch(
            "dtocean_electrical.inputs.ElectricalArrayData._set_footprints",
            return_value=None)
            
    array_output = [0.2, 0.8]
    
    result = ElectricalArrayData(machine_data,
                                 None,
                                 None,
                                 1,
                                 array_output)
        
    assert result.machine_data.power_factor == [(0.25, 1.0), (0.75, 1.0)]


def test_ElectricalArrayData_set_histogram_edges_list(mocker):
    
    # Fake a ElectricalMachineData object
    machine_data = mocker.Mock()
    machine_data.power = 1e6
    machine_data.power_factor = [[0.5, 0],
                                 [0.8, 1],
                                 [1, 1]]
    
    mocker.patch(
            "dtocean_electrical.inputs.ElectricalArrayData._set_footprints",
            return_value=None)
            
    array_output = [0.2, 0.2, 0.2, 0.2, 0.2]
    
    result = ElectricalArrayData(machine_data,
                                 None,
                                 None,
                                 1,
                                 array_output)
    
    zeros = [x[1] for x in result.machine_data.power_factor[:3]]
    ones = [x[1] for x in result.machine_data.power_factor[3:]]
    
    assert set(zeros) == set([0])
    assert set(ones) == set([1])
