#!/usr/bin/env ruby

require 'open-uri'

('GSM2139471'..'GSM2139662').each do |gsm|
  response = open("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=#{gsm}")
  disease_match = response.read.match /disease: ([A-Z]*)/
  puts "#{gsm}\t#{disease_match[1]}"
end