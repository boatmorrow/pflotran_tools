'''
make a pretty picture of the source concentrations
'''

from pflotran_tools import *

dd = ReadPflotranObsPt(obsptfile='observation-0.tec');

plot(dd['Time [y]'],dd['Total Am_241(aq)'],label='Am241');
plot(dd['Time [y]'],dd['Total Np_237(aq)'],label='Np237');
semilogy(dd['Time [y]'],dd['Total U_233(aq)'],label='U233');
semilogy(dd['Time [y]'],dd['Total Th_229(aq)'],label='Th229');
semilogy(dd['Time [y]'],dd['Total I_129(aq)'],label='I129');

xlabel('Conc [M]');
ylabel('Time (yr)');
legend(loc='best');
show()
