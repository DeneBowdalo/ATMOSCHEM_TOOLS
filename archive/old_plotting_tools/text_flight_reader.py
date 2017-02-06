import csv
import numpy as np
import datetime

def reader_text(filename):
    year_month_day_hours_mins_2 = []
    processed_time=[]
    update = 0
    for files in filename:
        print files
        reader=csv.reader(open(files,'rU'), delimiter='\t')
        countnum =0
        for row in reader:
            if row[0] != 'time0':
                new = row[1:] 
                try:
                    data.append(new)
                except:
                    data=[new]

            if row[0] == 'time0':
                names = row[1:4]
                raw_date = row[5]
                year = raw_date[:4]
                month = raw_date[5:7]
                day = raw_date[8:10]
                hour = raw_date[11:13]
                mins = raw_date[14:16]
                secs = raw_date[17:19]
            countnum+=1 
        base_time= year+month+day+hour+mins+secs
        processed_time.append(base_time)
        secs = int(secs)
        mins = int(mins)
        hour = int(hour)
        day  = int(day)
        month= int(month)
        year = int(year)
        for i in range(1,countnum,1):  
            if i > 1:

                a = datetime.datetime(year,month,day,hour,mins,secs)
            
                b = a + datetime.timedelta(0,1)
                c = str(a)
                d = str(b)
                year_p = d[:4]
                month_p =d[5:7]
                day_p = d[8:10]
                hour_p = d[11:13]
                mins_p = d[14:16]
                secs_p = d[17:]
                time_p = year_p+month_p+day_p+hour_p+mins_p+secs_p 
                processed_time.append(time_p)
                secs = int(secs_p)
                mins = int(mins_p)
                hour = int(hour_p)
                day  = int(day_p)
                month= int(month_p)
                year = int(year_p)
    
            else:
                a = datetime.datetime(year,month,day,hour,mins,secs)
  
    processed_time=np.array(processed_time) 
   
    
    for row in data:
        for ix, char in enumerate(row):
            if char == '':
                row[ix] = '000.00'
    data=np.array(data)
    data=np.float64(data)
   
    data= np.column_stack((data,processed_time))

    return data, names



