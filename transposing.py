import csv


#a = zip(*csv.reader(open("/home/james/summer2021/THADSNBit/30alpha.csv", "r")))
#csv.writer(open("ED30alphaT.csv",'w', newline='')).writerows(a)

#a = zip(*csv.reader(open("/home/james/summer2021/THADSNBit/40alpha.csv", "r")))
#csv.writer(open("ED40alphaT.csv",'w', newline='')).writerows(a)

a = zip(*csv.reader(open("/home/james/summer2021/THADSNBit/ED30alphaT.csv", "r")))
csv.writer(open("ED30alphaTT.csv",'w', newline='')).writerows(a)

a = zip(*csv.reader(open("/home/james/summer2021/THADSNBit/ED40alphaT.csv", "r")))
csv.writer(open("ED40alphaTT.csv",'w', newline='')).writerows(a)