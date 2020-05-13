#!/usr/bin/env ruby

r1_filename = ARGV[0] 
r2_filename = ARGV[1]                                                            
current_pair = 1

r1_lines = []
r2_lines = []

$stdin.each_line do |line|
  if line.start_with?("@HWI")
  	if line.chomp.end_with?("/1")
  	  current_pair = 1
    else
      current_pair = 2
    end
  end

  if current_pair == 1
    r1_lines.push(line)
  else
    r2_lines.push(line)
  end
end

File.write(r1_filename, r1_lines.join)
File.write(r2_filename, r2_lines.join)
