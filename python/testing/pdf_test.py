

x=np.linspace(0,10,1000)
plt.clf()
sigma=2
G = 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-x**2/(2*sigma**2))
M = x**2/np.sqrt(2*np.pi*sigma**2)*np.exp(-x**2/(2*sigma**2))
if 1:
    plt.plot(x,M/M.max(),label='M')
for fact in [0.5,1,2,7]:
    B=np.sinh(fact*x)
    BT=x*G*B
    plt.plot(x,BT/BT.max(),label='%0.f'%fact)
if 0:
    plt.plot(x,BT/BT.max(),label='BT')
    plt.plot(x,M/M.max(),label='M')
if 0:
    plt.plot(x,M,label='M')
    plt.plot(x,G,label='G')
#plt.plot(x,B)
    plt.plot(x,G*B,label='G*B')
    plt.plot(x,BT,label='x*G*B')
plt.legend(loc=0)
#plt.yscale('log')
plt.savefig('plots_to_sort/dorp.png')
