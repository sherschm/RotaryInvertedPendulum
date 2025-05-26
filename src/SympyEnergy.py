import sympy as sp

# Define symbolic variables and functions
t = sp.symbols('t')
theta1 = sp.Function('theta1')(t)
theta2 = sp.Function('theta2')(t)
theta1d = sp.Derivative(theta1, t)
theta2d = sp.Derivative(theta2, t)

# Define shorthand for sin/cos and powers
s1 = sp.sin(theta1)
c1 = sp.cos(theta1)
s2 = sp.sin(theta2)
c2 = sp.cos(theta2)

s1_2 = s1**2
c1_2 = c1**2
s2_2 = s2**2
c2_2 = c2**2
s1_4 = s1**4
c1_4 = c1**4
s2_4 = s2**4
c2_4 = c2**4
c2_3 = c2**3


T=(3.3880467526614666e-5 * (s1_2)*(theta1d**2) +
    8.987784385681152e-5 * (c2_2)*(theta1d**2) + 
    2.1690092154358354e-5 * (theta2d**2)*(s2_2) +
    3.3880467526614666e-5 * (c1_2)*(theta1d**2) +
    5.423644324764609e-5 * (s1_2)*c2*theta2d*theta1d + 
    0.00015401840209960938 * (c2_3)*theta2d*theta1d +
    5.423644324764609e-5*c2*theta2d*(c1_2)*theta1d + 
    2.1694577299058437e-5 * (s1_2)*(c2_2)*(theta2d**2) +
    2.1694577299058437e-5 * (s1_2)*(theta1d**2)*(s2_2) + 
    0.00010269582271575928*(c2_4)*(theta2d**2) +
    2.1694577299058437e-5 * (c2_2)*(theta2d**2)*(c1_2) + 
    2.1694577299058437e-5 * (c1_2)*(theta1d**2)*(s2_2) + 
    0.00015401840209960938 * (s1_2)*c2*theta2d*theta1d*(s2_2) + 
    0.00015401840209960938*c2*theta2d*(c1_2)*theta1d*(s2_2) + 
    0.00019252300262451172*(s1_4)*(theta1d**2)*(s2_2) +
    0.00020539164543151856*(s1_2)*(c2_2)*(theta2d**2)*(s2_2) +
    0.00038504600524902344*(s1_2)*(c1_2)*(theta1d**2)*(s2_2) + 
    0.00020539164543151856*(c2_2)*(theta2d**2)*(c1_2)*(s2_2) 
    + 0.00019252300262451172*(c1_4)*(theta1d**2)*(s2_2) + 
    0.00010269582271575928*(s1_4)*(theta2d**2)*(s2_4) + 
    0.00020539164543151856*(s1_2)*(theta2d**2)*(c1_2)*(s2_4) + 
    0.00010269582271575928*(theta2d**2)*(c1_4)*(s2_4)
    )


# Simplify the expression
simplified_T = sp.simplify(T)

# Print the result
print("Simplified Expression:")
sp.pprint(simplified_T, use_unicode=True)

with open("simplified_expression_pretty.txt", "w", encoding="utf-8") as f:
    f.write(sp.pretty(simplified_T, use_unicode=True))

    