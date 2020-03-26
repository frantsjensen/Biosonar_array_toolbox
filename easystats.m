function [results]=easystats(x)

A=mean(x);
B=std(x);
C=min(x);
D=max(x);

results = [A,B,C,D];

%disp(['Mean: ' num2str(A) '+/-' num2str(B) ' (1 std)'])
%disp(['Median: ' num2str(median(x))])