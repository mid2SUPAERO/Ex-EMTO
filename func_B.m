% GEOMETRIC MATRIX B USED IN STIFFNESS MATRIX CALCULATION
function B = func_B()
    syms s t;
    N1 = (1-s)*(1-t)/4; N2 = (1+s)*(1-t)/4;
    N3 = (1+s)*(1+t)/4; N4 = (1-s)*(1+t)/4;
    Ns = [N1, N2, N3, N4];
    Bs = sym(zeros(3,2,4));
    for i = 1:4
        Bs(:,:,i) = [[diff(Ns(i),s) 0];[0 diff(Ns(i),t)];[diff(Ns(i),t) diff(Ns(i),s)]];
    end
    B = [Bs(:,:,1),Bs(:,:,2),Bs(:,:,3),Bs(:,:,4)];
    B = matlabFunction(B);
end