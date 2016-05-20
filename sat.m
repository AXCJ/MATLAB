function [ output_args ] = sat( input_args )

if abs(input_args) > 0.9
    output_args = sign(input_args);
else
    output_args = input_args/3;
end

