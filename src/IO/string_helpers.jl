function searchforexactheaderstring(header_string, file_strings)
    # find the next "loop_"
    ind = findfirst(xx->header_string==xx, file_strings)
    if typeof(ind) == Int # must not be empty.
        return ind
    end

    return 0
end

function searchforheaderstring(header_string::String, file_strings)
    # find the next "loop_"
    ind = findfirst(xx->occursin(header_string, xx), file_strings)
    if typeof(ind) == Int # must not be empty.
        return ind
    end

    return 0
end

function searchforheaderstring(header_strings::Vector{String}, file_strings)
    # find the next "loop_"
    ind = findfirst(xx->any( occursin(header_strings[i], xx) for i = 1:length(header_strings) ), file_strings)
    if typeof(ind) == Int # must not be empty.
        return ind
    end

    return 0
end

function findblock(file_strings,
                    header_phrase,
                    termination_phrase,
                    file_label::String)
    #
    ind = searchforheaderstring(header_phrase, file_strings)
    if !(ind > 0)
        msg_string = "Error processing $(header_phrase) in $(file_label)."
        println(msg_string)
        return Vector{String}(undef, 0)
    end

    text_buffer = file_strings[ind+1:end]

    # find block termination.
    stop_ind = searchforheaderstring(termination_phrase, text_buffer)
    if !(stop_ind > 0)
        msg_string = "Error processing $(termination_phrase) in $(file_label)."
        println(msg_string)
        return Vector{String}(undef, 0)
    end

    return text_buffer[1:stop_ind-1]
end
# config_path = "/home/roy/MEGAsync/inputs/NMR/configs/cs_config.txt"
# file_strings = readlines(config_path)
# header_phrase = "(-)-Norepinephrine"
# termination_phrase = "end"
# Q = findblock(file_strings, header_phrase, termination_phrase, config_path)
