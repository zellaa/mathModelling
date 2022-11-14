function [y_out] = solve_lorentz(y0,aleph,ba,ra,ta)
%Author: Zella Baig, Date: 30/12/20
%This function solves solves the lorentz equations using the odsolver
%function
%inputs:
%aleph, ba, ra: constants corresponding to a, b, r
%y0: the 3x1 array of initial y values
%ta: the vector containing all time values
%output: this outputs an Nx3 array, with the columns corresponding to y_i
%and the rows to their value at a given time
%
%Usage: y = solve_lorentz(y_init,a,b,r,t)

%initialise constants, time, and initial y values
a=aleph;  
r=ra;
b=ba;
t=ta;
y_init = y0;

%create the anonymous functions
f1=@(t,p,q,rr) a*(q-p);  
f2=@(t,p,q,rr) p*r-q-p.*rr;
f3=@(t,p,q,rr) p.*q-b*rr;

%solve the equations using these functions
[y1_out,y2_out,y3_out] = odsolver(t,y_init,f1,f2,f3);

%store the outputs in a single variable
y_out(:,1) = y1_out;
y_out(:,2) = y2_out;
y_out(:,3) = y3_out;
end

