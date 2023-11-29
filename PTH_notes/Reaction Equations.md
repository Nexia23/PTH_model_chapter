In model calculated in the beginning for Steady State consideration of erythropoisis
$LDH_{Blood}=\frac{[LDH] * \frac{ln(2)} {t_{hlf,LDH,d}} * Vol_{blood}}{[E]\cdot \frac{2*ln(2)}{t_{E,d}} +\underbrace{J_{R,d}}_{=0}}$ with $t_{hl,LDH,d}$ = 4      [3-5 Tage](https://www.medicoconsult.de/ldh/)
$k_{P,b}= P(s_{P,d}\cdot[Hb]+k_{R,d0})-P(\frac{2\cdot ln(2)}{t_{P,a}})$ 
#TODO Where are these values coming from?
$s_{P,d} = 0.00071535$   slope of P death increase
$k_{P,d0} =  0.48924947$ default death rate of Precusors

$\frac{dP}{dt}=k_{P,b}-P(s_{P,d}\cdot[Hb]+k_{R,d0})-P(\frac{2\cdot ln(2)}{t_{P,a}})$ 
$t_{P,a}=t_{P,a,0}-\frac{t_{R,a,max}}{1+e^{-k_{R,a}\cdot(Hkt-Hkt_0)}}$ with $t_{P,a,0}=11$

$\frac{dR}{dt}=P\cdot(\frac{2\cdot ln(2)}{t_{P,a}})-R\cdot \frac{2\cdot ln(2)}{t_{R,a}}-R\cdot M\cdot k_{inf}-R\cdot(k_{R,d}+s_{BH}\cdot J_{oiE}),k_{R,d}=0$
$t_{R,a}=t_{R,a,0}+\frac{t_{R,a,max}}{1+e^{-k_{R,a}\cdot(Hkt-Hkt_0)}}, t_{R,a0}=0$

$\frac{dE}{dt}=R\cdot k_{R,a}-E\cdot M\cdot k_{inf}-E\cdot(k_{E,d}+s_{BH}\cdot J_{oiE})$
$BH=s_{BH}\cdot J_{oiE,d}$  with $s_{BH}=k_{E,d}\cdot x$, with x as Bystander Haemolysis factor
Not implemented Idea Hill Kinetics for BH as this can reach BH=0, as the linear function above, which logistic function cannot:
$BH=\frac{J_{oiE,d}^h}{BH_{max}\cdot J_{oiE,d}^h+K_m^h}$
$\frac{diE}{dt}= M\cdot(E+R)-iE\cdot (k_{iE,d,0} + (1-k_{iE,pit,frac}) * \frac{k_{iE, ART_{max}}*(ART^{h_{ART}}}{ART^{h_{ART}}+ID_{50}^{h_{ART}}})-iE\cdot k_{iE,rupture}$
$k_{iE,d}=k_{iE,d,0} + (1-k_{iE,pit,frac}) * \frac{k_{iE, ART_{max}}*(ART^{h_{ART}}}{ART^{h_{ART}}+ID_{50}^{h_{ART}}}$

$\frac{dM}{dt}=16\cdot iE\cdot r_{iE,rupture}-M\cdot k_{M,d}$


$\frac{doiE}{dt}=iE\cdot k_{iE,p}-\sum\limits_{i} \frac{doiE_i}{dt}=>\frac{doiE}{dt}=oiE\cdot k_{oiE,d}$
$k_{iE,pit}=k_{iE,pit,0} + k_{iE,pit,frac} * \frac{k_{iE, ART_{max}}*ART^{h_{ART}}}{ART^{h_{ART}}+ID_{50}^{h_{ART}}}$

For linear chain trick 12 species of oiE and $k_{oiE,d}\approx 0.963$ both values from fitting gamma distribution  
Event at t=0 ART set to 400, as medication given to patient
$\frac{dART}{dt}=ART\cdot\frac{ln(2)}{t_{ART,decay}}$
#TODO: Flux for cell deaths which release LDH
$\frac{dLDH}{dt}=\frac{LDH_{Blood}}{V_{Blood}}\cdot \sum J_{d, of RBC cells}-LDH\cdot \frac{ln(2)}{t_{hl,LDH,d}}$
$LDH_{Blood}=140-280 \frac{U}{L}$ only value nothing on degradation [link](https://www.ncbi.nlm.nih.gov/books/NBK557536/?report=printable)


