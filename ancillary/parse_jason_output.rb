#!/usr/bin/env ruby

require 'json'

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
