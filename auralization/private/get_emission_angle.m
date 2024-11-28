function output = get_emission_angle(input)

for i = 1:size(input,1)

    output.phix(i) = input{i}.phix; 
    output.phiy(i) = input{i}.phiy; 

end

end

