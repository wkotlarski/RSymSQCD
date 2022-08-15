#!/usr/bin/env ruby

require 'json'

if ARGV.length != 1
  puts "Expected one argument (JSON file name). Received #{ARGV.length}"
  exit 1
end

if !File.file?(ARGV[0])
  puts "No such file #{ARGV[0]}"
  exit 1
end

content = JSON.parse(File.read(ARGV[0]))

xsec = 0.0
err = 0.0
content["cross sections"].each do |key, _|
  content["cross sections"][key].each do |key2, _|
    xsec += content["cross sections"][key][key2]["res"]
    err  += content["cross sections"][key][key2]["err"]**2
  end
end

err = Math.sqrt(err)

# error exponent in scientific notation
# e.g. 2.0e-04 => -4

err_as_string = '%.1e' % err
err_exp = err_as_string.to_s.split('e').last.to_i
xsec_as_string = '%.1e' % xsec
xsec_exp = xsec_as_string.to_s.split('e').last.to_i

zeros = Array.new(xsec_exp-err_exp, "0").insert(1, '.').join('')
err_as_str2 = (err*10**(xsec_exp-err_exp)).round.to_s

sign_of_exp = xsec_exp>=0 ? '+' : '-'
puts "#{xsec} +/- #{err} (#{sprintf('%.1e', 1e+2*err/xsec)} %)"
puts "%.#{err_exp.abs+1}e +/- #{zeros}#{err_as_str2}e#{sign_of_exp}0#{xsec_exp}" % xsec
