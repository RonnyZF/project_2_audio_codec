import csv
import matplotlib.pyplot as plt

d_input=[]
d_out=[]

#------------------------------------------------

with open('/var/nfs/rootfs/home/root/Input.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count=0
    for row in csv_reader:
        d=float(row[0])
        d_input.append(d)

with open('/var/nfs/rootfs/home/root/Decoder.csv') as csv_file:
#with open('/home/project2/eclipse-workspace/prueba_decoder_ubuntu/Decoder.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count=0
    for row in csv_reader:
        r=float(row[0])
        d_out.append(r)

print(len(d_out))

time = [i for i in range(len(d_input))]

fig=plt.figure()
ax = fig.add_subplot(111)
ax.plot(time, d_input, 'r', label='Input Signal', Lw=1)
ax.plot(time, d_out[0:len(time)], 'b', label='Output Signal', Lw=1)
ax.axis([10000,10600,-1.25,1.25])
ax.set_xlabel(r'time (s)')
ax.set_ylabel(r'Ampl (V)')
ax.grid(True)
#plt.tight_layout()
plt.title(r'Input and Output signals')
plt.legend(loc='lower right')
plt.savefig("señales.png")
plt.show()
