import matplotlib.pyplot as plt


labels = ['$96^3$', '$128^3$', '$160^3$', '$192^3$']
mem = [0.374, 0.802, 1.4, 2.4]
time = [48.7, 171.8, 293.9, 514.9]

plt.figure(figsize=(8, 4))
plt.subplot(121)
bars=plt.bar(labels, time, color='gray', label=labels)
plt.bar_label(bars)
plt.xlabel('Grid size')
plt.ylabel('Runtime / seconds')
plt.title('(a)')


plt.subplot(122)


bars=plt.bar(labels, mem,  color='gray', label=labels)
plt.bar_label(bars)
plt.xlabel('Grid size')
plt.ylabel('Memory / GB')
plt.title('(b)')

plt.tight_layout(pad=0.5)
plt.savefig('runtime_memory.png')
plt.show()
