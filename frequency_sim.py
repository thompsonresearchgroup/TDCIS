import math
import sys


def read_parameters_from_file(filename):
    with open(filename, 'r') as file:
        params = file.readline().strip().split()
        number_of_atoms = int(params[0])
        delta_time = float(params[1])
        vibration_frequency = float(params[2])
        vibration_amplitude = float(params[3])

        return number_of_atoms, delta_time, vibration_frequency, vibration_amplitude

def get_last_line(filename):
    """Reads the specified file and returns its last line."""
    with open(filename, 'r') as file:
        non_empty_lines = [line for line in file.readlines() if line.strip()]
        return non_empty_lines[-1] if non_empty_lines else ""

def delta_x_y_z(file_name, time, delta_t, vibration_frequency, amplitude, move_nuclear_params=True):
    w = vibration_frequency * 7.25163E-7 * 10000

    # Get the original last line if move_nuclear_params is True
    original_last_line = ""
    if move_nuclear_params:
        original_last_line = get_last_line(file_name)
        print("Original last line:", original_last_line)

    processed_lines = []

    with open('d_xyz.txt', 'r') as file_input:
        lines = file_input.readlines()
        for line in lines:
            if lines.index(line) < 2:
                processed_lines.append(line.strip())
                first_value = int(processed_lines[0].split()[0])
                processed_lines[0] = str(first_value)
                continue

            values = line.split()
            if len(values) < 4:
                processed_lines.append(line.strip())
                continue

            charge, x, y, z = values
            dx = amplitude * float(x) * math.sin(w * (time + delta_t)) - amplitude * float(x) * math.sin(w * time)
            dy = amplitude * float(y) * math.sin(w * (time + delta_t)) - amplitude * float(y) * math.sin(w * time)
            dz = amplitude * float(z) * math.sin(w * (time + delta_t)) - amplitude * float(z) * math.sin(w * time)

            output_line = "{0} {1} {2} {3}".format(charge, dx, dy, dz)
            processed_lines.append(output_line)

    if move_nuclear_params and original_last_line:
        processed_lines[-1] = original_last_line.strip()

    with open(file_name, 'w') as f_output:
        for line in processed_lines:
            f_output.write(line + "\n")


if __name__ == '__main__':
    if len(sys.argv) == 1 or sys.argv[1] == '--help':
        print('')
        print('|--------------------------------------------------------------------------------------------------------------------')
        print('| Usage: python script.py file time delta_t output_file vibration_frequency amplitude                                |')
        print('|--------------------------------------------------------------------------------------------------------------------')
        print('|                                                                   ')
        print('|  script.py - refers to the name of this python file')
        print('|')
        print('|  file - refers to the frequencies file containing displacement dx, dy, and dz and is labelled as originald_xyz.txt')
        print('|')
        print('|     - its format is:')
        print('|       3 ----> number of atoms')
        print('|       h2q ----> a string to label they type of atoms: 2 h atoms and a ghost atom q')
        print('|       0.0 0.0 0.0 0.0  ----> the first value is the charge, the subsequent values are dx, dy, dz')
        print('|       0.0 0.0 0.0 0.0  ----> values of dx, dy, and dz are acquired from .log gaussian freq calculation')
        print('|       0.0 0.0 0.0 0.0  ')
        print('|  time - time (t) during simulation')
        print('|  delta-t = time-step == nucleat-time')
        print('|  output_file = file to which new dimensions are written and which allows us move the point charge -> nuc.txt')
        print('|  vibration frequency - acquired from gaussian .log file after freq calculation.')
        print('|  vibration amplitude - define as necessary and realistic')
        print('|')
        print('|____________________________________________________________________________________________________________________')
        print('')
        sys.exit(0)
    elif len(sys.argv) != 3:
        print(len(sys.argv))
        print('Error: Invalid number of arguments.')
        print('Usage: python script.py file time delta_t output_file vibration_frequency amplitude [move_nuclear_params]')
        sys.exit(1)

    file_name = sys.argv[1]
    time = float(sys.argv[2])


    # Check if the optional argument 'move_nuclear_params' is present and true
    if len(sys.argv) == 8 and sys.argv[7].lower() == 'true':
        move_nuclear_params = True
    else:
        move_nuclear_params = False
    number_of_atoms, delta_t, vibration_frequency, vibration_amplitude = read_parameters_from_file("d_xyz.txt")
    delta_x_y_z(file_name, time, delta_t, vibration_frequency, vibration_amplitude, move_nuclear_params=True)
