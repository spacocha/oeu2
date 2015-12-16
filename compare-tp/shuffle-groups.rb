#!/usr/bin/env ruby
#
# author: scott w olesen <swo@mit.edu>

require 'set'
require 'arginine'

par = Arginine::parse do
  desc "shuffle groups and look for co-occurences"
  opt :n, help: "number of times to do each kind of shuffling", default: 10
  arg :groups1, help: "first groups file"
  arg :groups2, help: "second groups file"
  flag :verbose, help: "show information about keys?"
end

$dat1 = File.readlines(par[:groups1]).drop(1).map { |l| l.chomp.split }.to_h
$dat2 = File.readlines(par[:groups2]).drop(1).map { |l| l.chomp.split }.to_h

# find common keys
keys1 = Set.new $dat1.keys
keys2 = Set.new $dat2.keys
keys = keys1.intersection(keys2).to_a
xor_keys = keys1 ^ keys2

if par[:verbose]
  $stderr.puts "there are #{keys.length} common ids and #{xor_keys.length} xor ids"
end

# rework into simple arrays with the common keys
group1 = keys.map { |k| $dat1[k] }
group2 = keys.map { |k| $dat2[k] }

def pair_type(i, j)
  case [group1[i] == group1[j], group2[i] == group2[j]]
  when [true, true]
    :both
  when [true, false], [false, true]
    :one
  when [false, false]
    :neither
  end
end

def pair_types(g1, g2)
  raise unless g1.length == g2.length
  idx = (1..(g1.length)).to_a
  dat = {both: 0, one: 0, neither: 0}
  idx.combination(2).each do |i, j|
    case [g1[i] == g1[j], g2[i] == g2[j]]
    when [true, true]
      dat[:both] += 1
    when [true, false], [false, true]
      dat[:one] += 1
    when [false, false]
      dat[:neither] += 1
    end
  end
  dat
end

def display_pairs(type, g1, g2)
  x = pair_types(g1, g2)
  puts [type, x[:both], x[:one], x[:neither]].join("\t")
end

# show the original results
puts %w<type both one neither>.join("\t")
display_pairs("base", group1, group2)
n = 10

1.upto(n).each { |i| display_pairs("shuffle1", group1.shuffle, group2) }
1.upto(n).each { |i| display_pairs("shuffle2", group1, group2.shuffle) }
1.upto(n).each { |i| display_pairs("shuffle_both", group1.shuffle, group2.shuffle) }
