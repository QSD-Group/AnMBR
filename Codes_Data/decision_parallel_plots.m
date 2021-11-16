clear 'individual_vals','sum_vals';
i=1;
individual_vals = zeros(100000,6);
sum_vals = zeros(100000,6);
for a=1:length(A)
    for a1=1:length(A1)
        for b=1:length(B)
            for c=1:length(C)
                for d=1:length(D)
                    for h=1:length(H)
                        individual_vals(i,1)=A(a);
                        individual_vals(i,2)=A1(a1);
                        individual_vals(i,3)=B(b);
                        individual_vals(i,4)=C(c);
                        individual_vals(i,5)=D(d);
                        individual_vals(i,6)=H(h);
                        sum_vals(i,1)=A(a);
                        sum_vals(i,2)=A(a)+A1(a1);
                        sum_vals(i,3)=A(a)+A1(a1)+B(b);
                        sum_vals(i,4)=A(a)+A1(a1)+B(b)+C(c);
                        sum_vals(i,5)=A(a)+A1(a1)+B(b)+C(c)+D(d);
                        sum_vals(i,6)=A(a)+A1(a1)+B(b)+C(c)+D(d)+H(h);
                        i=i+1;
                    end
                end
            end
        end
    end
end
