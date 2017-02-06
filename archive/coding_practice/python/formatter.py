def format_2dp(number):
    print('Number is {:.2f}'.format(number))
   
import fileinput

for number in fileinput.input():
    float_num = float(number)
    format_2dp(float_num)
