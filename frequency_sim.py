import fileinput
import math
import sys

def delta_x_y_z(file_name, time, delta_t, output_file, vibration_frequency, amplitude):
    """
    1 second = 4.134137E16 atomic units
    1 cm^-1 = 2.997925E10 Hz = 2.997925E10 s^-1
    1 cm^-1 = 2.997925E10/4.134137E16 = 7.25163E10-7 a.u.^-1
    """
    w = vibration_frequency *  7.25163E10-7

    # Open the output file for writing
    with open(output_file, 'w') as f_output:

        # Iterate over the input lines using fileinput
        for line in fileinput.input(files=file_name):

            # Write the first two lines from the original file
            if fileinput.isfirstline():
                f_output.write(line)
                continue

            values = line.split()

            # Make sure the line has at least four values
            if len(values) < 4:
                f_output.write(line)
                continue

            charge, x, y, z = values

            dx = amplitude*float(x)*math.sin(w*(time + delta_t)) - amplitude*float(x)*math.sin(w*time)
            dy = amplitude*float(y)*math.sin(w*(time + delta_t)) - amplitude*float(y)*math.sin(w*time)
            dz = amplitude*float(z)*math.sin(w*(time + delta_t)) - amplitude*float(z)*math.sin(w*time)

            # Write the updated coordinates to the output file
            f_output.write(str(charge) + ' ' + str(dx) + ' ' + str(dy) + ' ' + str(dz) + '\n')


if __name__ == '__main__':
    if len(sys.argv) == 1 or sys.argv[1] == '--help':
        print('')
        print('--------------------------------------------------------------------------------------')
        print('| Usage: python script.py file time delta_t output_file vibration_frequency amplitude |')
        print('| Guide: python ~/TDCIS/frequency_sim.py - enclose this using single quotation marks |')
        print('--------------------------------------------------------------------------------------')
        print('')
        sys.exit(0)
    elif len(sys.argv) != 7:
        print('Error: Invalid number of arguments.')
        print('Usage: python script.py file time delta_t output_file vibration_frequency amplitude')
        sys.exit(1)

    file_name = sys.argv[1]
    time = float(sys.argv[2])
    delta_t = float(sys.argv[3])
    output_file = sys.argv[4]
    vibration_frequency = float(sys.argv[5])
    amplitude = float(sys.argv[6])

    # Call delta_x_y_z with the specified arguments
    delta_x_y_z(file_name, time, delta_t, output_file, vibration_frequency, amplitude)

