function NewtonInterpolation(data, c)

syms x;

dataX = data(:, 1);
dataY = data(:, 2);

% Create divided diffrence table
n = length(data);
f = zeros(n, n);    % Divided diffrence table
f(:, 1) = dataY;
for j = 1:n
    for i = 1:n - j
        f(i, j + 1) = (f(i + 1, j) - f(i, j))/(dataX(i + j) - dataX(i));
    end
end

% Find start index
delta = abs(c - dataX);
delta_min = min(delta);
delta_min_index = find(delta == delta_min);
start_index = delta_min_index(1);

% Using forward formula or backward formula???
delta1 = abs(start_index - [1 n]);

if delta1(1) <= delta1(2)
    % Newton’s forward interpolation
    P = f(start_index, 1);
    for i = 2:n - start_index + 1
        t = 1;
        for k = start_index:start_index + i - 2
            t = t*(x - dataX(k));
        end
        P = P + t*f(start_index, i);
    end
else
    % Newton’s backward interpolation
    P = f(start_index, 1);
    for i = 2:start_index
        t = 1;
        for k = start_index:-1:start_index - i + 2
            t = t*(x - dataX(k));
        end
        P = P + t*f(start_index - i + 1, i);
    end
end
P = simplify(P);
disp(P);
P_c = double(subs(P, c));
disp(P_c);
end



