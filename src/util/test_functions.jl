function writealpha()
    str = ""
    for key in sort(collect(keys(coordinates)))
        str *= key
    end

    writeSentence(str)    
end