function [x_predict, y_predict] = localise(stack, u, v)
    third = round(length(stack)/3);
    sublens = stack((v+1)*third+1 : (v+2)*third, (u+1)*third+1 : (u+2)*third);
    temp = sublens;
    
    surrounding = 4; % surrounding pixels
    x_predict = [];
    y_predict = [];

    maximum = max(sublens, [], 'all');
    [x_max, y_max] = find(sublens == maximum);
    % first point
    x_max_1 = x_max(1);
    y_max_1 = y_max(1);

    % second point
    if u == 0 && v == 0
        x_max_2 = x_max_1;
        y_max_2 = y_max_1;
        
    else
        if length(x_max) > 1 % if more than 1 maximum (2)
            x_max_2 = x_max(2);
            y_max_2 = y_max(2);
        else
            temp = sublens(sublens < maximum);
            maximum = max(temp, [], 'all');
            [x_max, y_max] = find(sublens == maximum);
            x_max_2 = x_max(1);
            y_max_2 = y_max(1);
        end

        while true
            if sqrt((x_max_1-x_max_2)^2 + (y_max_1-y_max_2)^2) < 2 % within 1 pixels
                temp = temp(temp < maximum);
                maximum = max(temp, [], 'all');
                [x_max, y_max] = find(sublens == maximum);
                x_max_2 = x_max(1);
                y_max_2 = y_max(1);
            else
                break;
            end
        end
    end

    subsection_1 = sublens(x_max_1-surrounding:x_max_1+surrounding, y_max_1-surrounding:y_max_1+surrounding);
    xx_1 = 1 : length(subsection_1);
    yy_1 = 1 : length(subsection_1);
    [fitresult, ~, ~, ~, ~, ~] = fmgaussfit(xx_1, yy_1, subsection_1);
    x_predict(1) = fitresult(5) + x_max_1 - surrounding - 1;
    y_predict(1) = fitresult(6) + y_max_1 - surrounding - 1;

    subsection_2 = sublens(x_max_2-surrounding:x_max_2+surrounding, y_max_2-surrounding:y_max_2+surrounding);
    xx_2 = 1 : length(subsection_2);
    yy_2 = 1 : length(subsection_2);
    [fitresult, ~, ~, ~, ~, ~] = fmgaussfit(xx_2, yy_2, subsection_2);
    x_predict(2) = fitresult(5) + x_max_2 - surrounding - 1;
    y_predict(2) = fitresult(6) + y_max_2 - surrounding - 1;

    temp = x_predict;
    x_predict = y_predict;
    y_predict = temp;

    x_predict = x_predict + (u + 1) * third; % sublens -> stack
    y_predict = y_predict + (v + 1) * third;
    
    x_predict = x_predict';
    y_predict = y_predict';
    
end

