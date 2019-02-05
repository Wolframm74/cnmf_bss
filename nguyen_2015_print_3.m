function nguyen_2015_print_3(input_arg)

n_current=input_arg(1);
L_value_3rd_combined_costfun=input_arg(3);

    if (input_arg(2)==1)

        disp('Just updated Xhat, B, and C.')
        display(n_current)
        display(L_value_3rd_combined_costfun)

    elseif (input_arg(2)==2)

        disp('Just updated T.')
        display(n_current)
        %display(L_value)

    elseif (input_arg(2)==3)

        disp('Just updated V.')
        display(n_current)
        %display(L_value)

    elseif (input_arg(2)==4)

        disp('Just updated W.')
        display(n_current)
        %display(L_value)

    elseif (input_arg(2)==5)

        disp('Just updated Phi_S.')
        display(n_current)
        %display(L_value)

    elseif (input_arg(2)==6)

        disp('Just updated Z.')
        display(n_current)
        %display(L_value)
        
    elseif (input_arg(2)==7)

        disp('Just updated Y.')
        display(n_current)
        %display(L_value)        
        
    end

%pause


end
