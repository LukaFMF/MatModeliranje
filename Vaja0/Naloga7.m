1;
function d = odvod(p)
    dolz = size(p)(2);
    stopnja = dolz - 1;
    d = zeros(1,stopnja);
    for i = 1:dolz - 1
        d(i) = p(i)*stopnja;
        stopnja--;
    end
end

odvod([3 1 1 1])