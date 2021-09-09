
# for i in range(100):
#    print()
#    print('        <volume name ="ShieldVol_' + str(i) + '">')
#    print('            <materialref ref="' + Materials[i] + '"/>')
#    print('            <solidref ref="Shield_' + str(i) + '"/>')
#    print('        </volume>')

# for y in range(32):
#     for x in range(32):
#         print(y*x)
#         print("/detector/add Sivol_" + str(x+y*32))
#         print("A:", (31 - x) * (31 - y), "B:", x * (31 - y), "C:", (31 - x) * y, "D:", x * y)
#         print("Sum", (31 - x) * (31 - y) + x * (31 - y) + (31 - x) * y + x * y)

# x = 10
# y = 0
# print("/detector/add Sivol_" + str(x+y*32))
# print("A:", (31 - x) * (31 - y), "B:", x * (31 - y), "C:", (31 - x) * y, "D:", x * y)
# print("Sum", (31 - x) * (31 - y) + x * (31 - y) + (31 - x) * y + x * y)



#for i in range(1296):
#    print("/detector/add Sivol_" + str(i))

DensityLead = 11.35
DensityTungsten = 19.3
DensityHDPE = 0.94

for i in range(9):
    print('<material name="PE-Pb-' + str((i+1)*10) + '">')
    print('     <D value="' + str(DensityHDPE*(9-i)/10 + DensityLead*(i+1)/10) + '" unit="g/cm3"/>')
    print('     <fraction n="' + str((9-i)/10) + '" ref="Polyethylene"/>')
    print('     <fraction n="' + str((i+1)/10) + '" ref="lead"/>')
    print('</material>')


