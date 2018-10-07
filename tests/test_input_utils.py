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

import pytest
import pandas as pd

from dtocean_electrical.input_utils.utils import set_burial_from_bpi
from dtocean_electrical.input_utils.input_tests import check_inputs


def test_check_inputs_error_technology(mocker):
    
    # Fake a database
    df_dict = {"id": [1, 2],
               "v_rate": [6600, 11000],
               "a_air": [50, 50]}
    
    df = pd.DataFrame(df_dict)

    input_data = mocker.Mock()
    input_data.array_data.machine_data.technology = "flying"
    input_data.array_data.machine_data.power = 200000
    input_data.array_data.machine_data.connection = 'dry-mate'
    input_data.array_data.machine_data.voltage = 6600
    input_data.array_data.n_devices = 10
    input_data.options.user_umbilical = 1
    input_data.database.dynamic_cable = df
    input_data.array_data.machine_data.max_current =  20
    
    errstatus, errstr = check_inputs(input_data)
    
    assert errstatus == -1
    assert len(errstr) == 1
    
    
def test_check_inputs_error_power(mocker):
    
    # Fake a database
    df_dict = {"id": [1, 2],
               "v_rate": [6600, 11000],
               "a_air": [50, 50]}
    
    df = pd.DataFrame(df_dict)

    input_data = mocker.Mock()
    input_data.array_data.machine_data.technology = "fixed"
    input_data.array_data.machine_data.power = 1
    input_data.array_data.machine_data.connection = 'dry-mate'
    input_data.array_data.machine_data.voltage = 6600
    input_data.array_data.n_devices = 10
    input_data.options.user_umbilical = 1
    input_data.database.dynamic_cable = df
    input_data.array_data.machine_data.max_current =  20
    
    errstatus, errstr = check_inputs(input_data)
    
    assert errstatus == -1
    assert len(errstr) == 1


def test_check_inputs_error_connection(mocker):
    
    # Fake a database
    df_dict = {"id": [1, 2],
               "v_rate": [6600, 11000],
               "a_air": [50, 50]}
    
    df = pd.DataFrame(df_dict)

    input_data = mocker.Mock()
    input_data.array_data.machine_data.technology = "fixed"
    input_data.array_data.machine_data.power = 200000
    input_data.array_data.machine_data.connection = 'best-mate'
    input_data.array_data.machine_data.voltage = 6600
    input_data.array_data.n_devices = 10
    input_data.options.user_umbilical = 1
    input_data.database.dynamic_cable = df
    input_data.array_data.machine_data.max_current =  20
    
    errstatus, errstr = check_inputs(input_data)
    
    assert errstatus == -1
    assert len(errstr) == 1
    
    
def test_check_inputs_error_voltage(mocker):
    
    # Fake a database
    df_dict = {"id": [1, 2],
               "v_rate": [6600, 11000],
               "a_air": [50, 50]}
    
    df = pd.DataFrame(df_dict)

    input_data = mocker.Mock()
    input_data.array_data.machine_data.technology = "floating"
    input_data.array_data.machine_data.power = 200000
    input_data.array_data.machine_data.connection = 'dry-mate'
    input_data.array_data.machine_data.voltage = 3300
    input_data.array_data.n_devices = 10
    input_data.options.user_umbilical = 1
    input_data.database.dynamic_cable = df
    input_data.array_data.machine_data.max_current =  20
    
    errstatus, errstr = check_inputs(input_data)
    
    assert errstatus == -1
    assert len(errstr) == 1


def test_check_inputs_error_arraypower(mocker):
    
    # Fake a database
    df_dict = {"id": [1, 2],
               "v_rate": [6600, 11000],
               "a_air": [50, 50]}
    
    df = pd.DataFrame(df_dict)

    input_data = mocker.Mock()
    input_data.array_data.machine_data.technology = "floating"
    input_data.array_data.machine_data.power = 200000
    input_data.array_data.machine_data.connection = 'dry-mate'
    input_data.array_data.machine_data.voltage = 6600
    input_data.array_data.n_devices = -1
    input_data.options.user_umbilical = 1
    input_data.database.dynamic_cable = df
    input_data.array_data.machine_data.max_current =  20
    
    errstatus, errstr = check_inputs(input_data)
    
    assert errstatus == -1
    assert len(errstr) == 1
    
    
def test_check_inputs_error_max_current(mocker):
    
    # Fake a database
    df_dict = {"id": [1, 2],
               "v_rate": [6600, 11000],
               "a_air": [50, 50]}
    
    df = pd.DataFrame(df_dict)

    input_data = mocker.Mock()
    input_data.array_data.machine_data.technology = "floating"
    input_data.array_data.machine_data.power = 200000
    input_data.array_data.machine_data.connection = 'dry-mate'
    input_data.array_data.machine_data.voltage = 6600
    input_data.array_data.n_devices = 10
    input_data.options.user_umbilical = 1
    input_data.database.dynamic_cable = df
    input_data.array_data.machine_data.max_current =  100
    
    errstatus, errstr = check_inputs(input_data)
    
    assert errstatus == -1
    assert len(errstr) == 1


def test_check_inputs(mocker):
    
    # Fake a database
    df_dict = {"id": [1, 2],
               "v_rate": [6600, 11000],
               "a_air": [50, 50]}
    
    df = pd.DataFrame(df_dict)

    input_data = mocker.Mock()
    input_data.array_data.machine_data.technology = "floating"
    input_data.array_data.machine_data.power = 200000
    input_data.array_data.machine_data.connection = 'dry-mate'
    input_data.array_data.machine_data.voltage = 6600
    input_data.array_data.n_devices = 10
    input_data.options.user_umbilical = 1
    input_data.database.dynamic_cable = df
    input_data.array_data.machine_data.max_current =  20
    
    errstatus, errstr = check_inputs(input_data)
    
    assert errstatus == 0
    assert len(errstr) == 0


@pytest.mark.parametrize("test_input, expected", [
    ('hard glacial till', 0),
    ('cemented', 0),
    ('soft rock coral', 0),
    ('hard rock', 0),
    ('gravel cobble', 0),
    ('loose sand', 0.5),
    ('medium sand', 0.5),
    ('dense sand', 0.5),
    ('very soft clay', 1.0),
    ('soft clay', 1.0),
    ('firm clay', 1.0),
    ('stiff clay', 1.0)
    ])
def test_set_burial_from_bpi(test_input, expected):
    
    row = {'layer 1 type': test_input}
    bpi = set_burial_from_bpi(row)
    
    assert bpi == expected


def test_set_burial_from_bpi_error():
    
    row = {'layer 1 type': "shingle"}
    
    with pytest.raises(ValueError):
        set_burial_from_bpi(row)
