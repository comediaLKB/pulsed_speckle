function output = circ_mask(input_size,center,radius)

output = zeros(input_size);

for idx1 = 1:input_size(1)
    for idx2 = 1:input_size(2)

        if sqrt((idx1-center(1)).^2+(idx2-center(2)).^2) <= radius
            output(idx1,idx2) = 1;
        end
        
    end
end

end