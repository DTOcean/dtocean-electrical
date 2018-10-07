# -*- coding: utf-8 -*-

#    Copyright (C) 2017-2018 Mathew Topper
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
