from netCDF4 import Dataset
import re

rootgrp = Dataset('/glade/work/edliu/models/CMAQ_v5.3.2/data/output_CCTM_v532_intel_Bench_2016_12SE1/CCTM_CGRID_v532_intel_Bench_2016_12SE1_20160701.nc', "r")
attr_dict = rootgrp.__dict__  # show all the attributes in the cgrid file
varstr = attr_dict['VAR-LIST']  # select the variable list in the map, this is a string
varstr = re.sub("\s+", ",", varstr.strip())  # replace whitespaces with comma

lst = varstr.split(',')
print(f'The index of NO2 is {lst.index("NO2")}')
print(f'The index of N2O5 is {lst.index("N2O5")}')
print(f'The index of HNO3 is {lst.index("HNO3")}')
print(f'The index of HCL is {lst.index("HCL")}')
print(f'The index of ASO4J is {lst.index("ASO4J")}')
print(f'The index of ASO4I is {lst.index("ASO4I")}')
print(f'The index of ANH4J is {lst.index("ANH4J")}')
print(f'The index of ANH4I is {lst.index("ANH4I")}')
print(f'The index of ANO3J is {lst.index("ANO3J")}')
print(f'The index of ANO3I is {lst.index("ANO3I")}')
print(f'The index of AECJ is {lst.index("AECJ")}')
print(f'The index of AECI is {lst.index("AECI")}')
print(f'The index of AH2OJ is {lst.index("AH2OJ")}')
print(f'The index of AH2OI is {lst.index("AH2OI")}')
print(f'The index of ANAJ is {lst.index("ANAJ")}')
print(f'The index of ACLJ is {lst.index("ACLJ")}')
print(f'The index of ACLI is {lst.index("ACLI")}')
print(f'The index of NH3 is {lst.index("NH3")}')
