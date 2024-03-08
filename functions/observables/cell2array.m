function output =  cell2array(input, row)

[~, col] = size( input ) ;

for i=1:row
    output = zeros(length(input{row}), col);
    for j=1:col
        output(:,j) = input{row, j};
    end
end


end