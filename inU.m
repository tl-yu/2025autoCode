function y = inU(P,u1,u2)
judge = zeros(3,3);
B1 = [-0.3 0.1 0.1;
    0.2 -0.2 -0.3;
    0.1 0.1 0.2];
B2 = [0.2 0.1 0.1;
    -0.3 -0.3 0.2;
    0.1 0.2 -0.3];
for i = 1:3
    for j = 1:3
        if (P(i,j)+u1*B1(j,i)+u2*B2(j,i) >= 0)&&(P(i,j)+u1*B1(j,i)+u2*B2(j,i) <= 1)
            judge(i,j) = 1;
        else
            judge(i,j) = 0;
        end
    end
end

if sum(sum(judge)) == 9
    y = 1;
else
    y = 0;
end

end