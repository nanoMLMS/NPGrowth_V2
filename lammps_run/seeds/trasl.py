
from snow.io.xyz import *
from snow.misc.rototranslation import *


el, coords = read_xyz("Cu1688_Cube.xyz")
print(len(el))

coords_trasl=translate_com_to_origin(coords, el)
write_xyz("Cu1688_Cube_com.xyz", el, coords_trasl)
