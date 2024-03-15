- A, B, C, D = symbols('A B C D', integer=True)
- T = symbols('T')
- p = solve($A*T^3 +B*T^2 + C*T + D$, T, dict=True)
- p is dict with three solutions
- print_latex(p[0])
$\left\{ T : - \frac{- \frac{3 C}{A} + \frac{B^{2}}{A^{2}}}{3 \sqrt[3]{\frac{\sqrt{- 4 \left(- \frac{3 C}{A} + \frac{B^{2}}{A^{2}}\right)^{3} + \left(\frac{27 D}{A} - \frac{9 B C}{A^{2}} + \frac{2 B^{3}}{A^{3}}\right)^{2}}}{2} + \frac{27 D}{2 A} - \frac{9 B C}{2 A^{2}} + \frac{B^{3}}{A^{3}}}} - \frac{\sqrt[3]{\frac{\sqrt{- 4 \left(- \frac{3 C}{A} + \frac{B^{2}}{A^{2}}\right)^{3} + \left(\frac{27 D}{A} - \frac{9 B C}{A^{2}} + \frac{2 B^{3}}{A^{3}}\right)^{2}}}{2} + \frac{27 D}{2 A} - \frac{9 B C}{2 A^{2}} + \frac{B^{3}}{A^{3}}}}{3} - \frac{B}{3 A}\right\}$
- repr(p[0])
python code: `{T: -(-3*C/A + B**2/A**2)/(3*(np.sqrt(-4*(-3*C/A + B**2/A**2)**3 + (27*D/A - 9*B*C/A**2 + 2*B**3/A**3)**2)/2 + 27*D/(2*A) - 9*B*C/(2*A**2) + B**3/A**3)**(1/3)) - (np.sqrt(-4*(-3*C/A + B**2/A**2)**3 + (27*D/A - 9*B*C/A**2 + 2*B**3/A**3)**2)/2 + 27*D/(2*A) - 9*B*C/(2*A**2) + B**3/A**3)**(1/3)/3 - B/(3*A)}`

- print_latex(p[1])
$\left\{ T : - \frac{- \frac{3 C}{A} + \frac{B^{2}}{A^{2}}}{3 \left(- \frac{1}{2} - \frac{\sqrt{3} i}{2}\right) \sqrt[3]{\frac{\sqrt{- 4 \left(- \frac{3 C}{A} + \frac{B^{2}}{A^{2}}\right)^{3} + \left(\frac{27 D}{A} - \frac{9 B C}{A^{2}} + \frac{2 B^{3}}{A^{3}}\right)^{2}}}{2} + \frac{27 D}{2 A} - \frac{9 B C}{2 A^{2}} + \frac{B^{3}}{A^{3}}}} - \frac{\left(- \frac{1}{2} - \frac{\sqrt{3} i}{2}\right) \sqrt[3]{\frac{\sqrt{- 4 \left(- \frac{3 C}{A} + \frac{B^{2}}{A^{2}}\right)^{3} + \left(\frac{27 D}{A} - \frac{9 B C}{A^{2}} + \frac{2 B^{3}}{A^{3}}\right)^{2}}}{2} + \frac{27 D}{2 A} - \frac{9 B C}{2 A^{2}} + \frac{B^{3}}{A^{3}}}}{3} - \frac{B}{3 A}\right\}$
- repr(p[1])
python code:`{T: -(-3*C/A + B**2/A**2)/(3*(-1/2 - sqrt(3)*I/2)*(sqrt(-4*(-3*C/A + B**2/A**2)**3 + (27*D/A - 9*B*C/A**2 + 2*B**3/A**3)**2)/2 + 27*D/(2*A) - 9*B*C/(2*A**2) + B**3/A**3)**(1/3)) - (-1/2 - sqrt(3)*I/2)*(sqrt(-4*(-3*C/A + B**2/A**2)**3 + (27*D/A - 9*B*C/A**2 + 2*B**3/A**3)**2)/2 + 27*D/(2*A) - 9*B*C/(2*A**2) + B**3/A**3)**(1/3)/3 - B/(3*A)}`
 
- print_latex(p[2])
$\left\{ T : - \frac{- \frac{3 C}{A} + \frac{B^{2}}{A^{2}}}{3 \left(- \frac{1}{2} + \frac{\sqrt{3} i}{2}\right) \sqrt[3]{\frac{\sqrt{- 4 \left(- \frac{3 C}{A} + \frac{B^{2}}{A^{2}}\right)^{3} + \left(\frac{27 D}{A} - \frac{9 B C}{A^{2}} + \frac{2 B^{3}}{A^{3}}\right)^{2}}}{2} + \frac{27 D}{2 A} - \frac{9 B C}{2 A^{2}} + \frac{B^{3}}{A^{3}}}} - \frac{\left(- \frac{1}{2} + \frac{\sqrt{3} i}{2}\right) \sqrt[3]{\frac{\sqrt{- 4 \left(- \frac{3 C}{A} + \frac{B^{2}}{A^{2}}\right)^{3} + \left(\frac{27 D}{A} - \frac{9 B C}{A^{2}} + \frac{2 B^{3}}{A^{3}}\right)^{2}}}{2} + \frac{27 D}{2 A} - \frac{9 B C}{2 A^{2}} + \frac{B^{3}}{A^{3}}}}{3} - \frac{B}{3 A}\right\}$
- repr(p[2])
python code:`{T: -(-3*C/A + B**2/A**2)/(3*(-1/2 + sqrt(3)*I/2)*(sqrt(-4*(-3*C/A + B**2/A**2)**3 + (27*D/A - 9*B*C/A**2 + 2*B**3/A**3)**2)/2 + 27*D/(2*A) - 9*B*C/(2*A**2) + B**3/A**3)**(1/3)) - (-1/2 + sqrt(3)*I/2)*(sqrt(-4*(-3*C/A + B**2/A**2)**3 + (27*D/A - 9*B*C/A**2 + 2*B**3/A**3)**2)/2 + 27*D/(2*A) - 9*B*C/(2*A**2) + B**3/A**3)**(1/3)/3 - B/(3*A)}`