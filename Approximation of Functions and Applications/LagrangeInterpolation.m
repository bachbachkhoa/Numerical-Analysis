function LagrangeInterpolation(dataList, c)
% Lagrange interpolation
% Input: 2-D array consists of 2 columns, the value column of x, 
%        and the value column of y corresponding to x; c in [x0, xn] 
% Output: polynomial interpolation P, P(c)

dataX = dataList(:, 1);
if c < dataX(1) | c > dataX(end)
    disp([num2str(x0), ' not in [', num2str(dataX(1)), ', ', num2str(dataX(end)), ']'])
else
    syms x;
    dataY = dataList(:, 2);
    w = prod(x - dataX);
    dw = diff(w);
    D = (x - dataX).*subs(dw, dataX);
    P = simplify(w*sum(dataY./D));
    P_c = double(subs(P, c));
    fprintf('\nP = ');
    disp(P);
    disp(['P(', num2str(c), ') = ', num2str(P_c)]);
end
end

% Example: dataList = [1 17; 2 17.5; 3 76; 4 210.5; 7 1970], c = 5
% P = 2*x^4 - 17*x^3 + 81*x^2 - (307*x)/2 + 209/2
% P(5) = 487
