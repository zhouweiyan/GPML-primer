x = [1;2;3;4;5];
y = x;
x_outliers = 3;
len = length(x);
for i = 1:length(x)
    i
    for j = 1:length(x_outliers)
        j
        if x(i)==x_outliers(j) 
%             len = len - 1;
            x(i) = []; y(i) = [];
        end
    end
    if i==length(x)
        break;
    end
    [x,y]
end