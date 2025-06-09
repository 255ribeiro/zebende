import platform
from pathlib import Path

def check_architecture():
    machine_type = platform.machine().lower()
    if machine_type == 'aarch64':
        return 'aarch64'
    elif machine_type in ['x86_64', 'amd64']:
        return "x86_64"
    else:
        return machine_type


def get_platform():

    sys_dict = {'windows': 'windows/zebendezig.dll', 'linux': 'linux/libzebendezig.so', 'darwin': 'macos/libzebendezig.dylib'}

    sys_info = platform.uname()

    lib_folder = Path(__file__).parent / 'zig_libs'

    lib_file = lib_folder / str(check_architecture() + "-" + sys_dict[sys_info.system.lower()])

    return lib_file


if __name__ == '__main__':
    print(get_platform())