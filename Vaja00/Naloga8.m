Delitelji(54)

function stD = Delitelji(m)
    if m > 1
        meja = floor(sqrt(m));
        stD = 2;
        for i = 2:meja
            if(mod(m,i) == 0)
                stD = stD + 2;
                if(i == meja)
                    stD = stD - 1;
                end
            end
        end
    else
        stD = 1;
    end
end