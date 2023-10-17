from GL import *

class timer():
    def __init__(self):
        pass

    def scrub_output_log(self,data_location):

            output_log = "%s/OutputLog"%data_location
            fptr = open(output_log)
            lines=fptr.readlines()
            fptr.close()

            cycle=[]
            time=[]
            frames=[]
            wall=[]

            for line in lines:
                lll = line.split()
                frame= int( lll[2][-4:])
                frames.append(frame)
                cycle.append(int(lll[3]))
                time.append(float(lll[4]))
                wall.append(float(lll[5]))

            cycle=nar(cycle)
            time=nar(time)
            wall=nar(wall)

            dt = (time[1:]-time[:-1])
            dc = (cycle[1:]-cycle[:-1])
            dwall=(wall[1:]-wall[:-1])
            dwall_dc = dwall/dc
            prob_dwall = dwall_dc+0
            prob_dwall.sort()
            mean_dwall = prob_dwall[:prob_dwall.size//2].mean()
            #mean_dwall = prob_dwall.mean()
            self.cycle=cycle
            self.time=time
            self.wall=wall
            self.dt=dt
            self.dc=dc
            self.dwall=dwall
            self.dtdc=dt/dc
            self.dwdc=dwall/dc
            self.rate_mask = np.abs(self.dwdc) < 10*mean_dwall

tim = timer()
tim.scrub_output_log("/data/cb1/Projects/P19_CoreSimulations/u202-Beta2")
plt.plot(tim.dtdc)
plt.savefig('plots_to_sort/dtdc')
plt.clf()
plt.plot(tim.dwdc[tim.rate_mask])
total_wall= tim.dwall[tim.rate_mask].sum()
plt.savefig('plots_to_sort/wall')

frame=118
ds = yt.load("/data/cb1/Projects/P19_CoreSimulations/u202-Beta2/GravPotential/DD%04d/data%04d"%(frame,frame))
ds.print_stats()
Nzones=45343488
seconds_per_update = tim.dwdc[-1]
zone_up_core_sec = Nzones*1/(8*seconds_per_update)
print(zone_up_core_sec)
