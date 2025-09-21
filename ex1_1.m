% made in matlab r2025a
% Exercise 1.1

% the calculation is based on the book page in the appendix of the
% assignment
lambda1to5 = [2.36602037;
               5.49780392;
               8.63937983;
              11.78097245;
              14.92256510];

lambda6plus = @(i) (4*i - 1)*(pi/4);

maxLambda = 6;
assert(maxLambda > 0,"Lambda must be larger than zero.")
if maxLambda <= 5
    lambda_i = lambda1to5(1:maxLambda,:);
else
    lambdasLeft = maxLambda - 5;
    lambda6toMax = zeros(lambdasLeft,1);
    for l = 1:lambdasLeft
        lambda6toMax(l,1) = lambda6plus(l+5);
    end
    lambda_i = [lambda1to5; lambda6toMax];

end