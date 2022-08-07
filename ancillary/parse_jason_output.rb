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

puts "#{xsec} +/- #{err}"
