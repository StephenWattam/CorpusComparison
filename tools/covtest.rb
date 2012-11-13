
corpus_one = {"common" => 1000,
              "normal" => 100,
              "rare" => 1,
              "weird2" => 23,
              "more" => 22,
              "and more" => 34 
}
corpus_two = {"common" => 100223,
              "normal" => 102,
              "rare" => 13,
              "weird" => 23
}

# 
# corpus_one = {
# 1 => 0.3716803,
# 2 => 0.2778111,
# 3 => 0.8152372,
# 4 => 0.7715097,
# 5 => 0.0163179,
# 6 => -0.4898738,
# 7 => -0.6060137,
# 8 => -0.8882970,
# 9 => 0.2913591,
# 10 => -0.3661791,
# 11 => 0.1320750,
# 12 => 0.2637229,
# 13 => -0.7390226,
# 14 => -0.0395929,
# 15 => 0.3387334,
# 16 => 0.8598541,
# 17 => 0.7388236,
# 18 => -0.5928083,
# 19 => 0.9226006,
# 10 => -0.3571427
# }
# 
# 
# corpus_two = {
# 1 => 0.6396969,
# 2 => 0.7942405,
# 3 => -0.6364473,
# 4 => -0.6845633,
# 5 => -0.6908862,
# 6 => -0.5034169,
# 7 => 0.5745298,
# 8 => -0.1247591,
# 9 => -0.5129564,
# 10 => 0.0745857,
# 11 => 0.0733665,
# 12 => -0.0118882,
# 13 => 0.1763471,
# 14 => 0.1027599,
# 15 => -0.9737805,
# 16 => 0.8747677,
# 17 => 0.9479392,
# 18 => 0.0843604,
# 19 => -0.3518961,
# 10 => -0.3034039
# }


stats_one = {}
stats_two = {}



# compute basic stats needed for other ops.
def stats(corpus, output)
  output[:mean] = 0.0
  output[:var]  = 0.0

  corpus.each{ |k,v|
    output[:mean] += v
  }
  output[:sum] = output[:mean]
  output[:mean] /= corpus.length

  corpus.each{ |k,v|
    output[:var] += (output[:mean] - v)**2
  }
  output[:var] /= corpus.length

  output[:sd] = output[:var].to_f ** 0.5

  return output
end


class CorpusComparison
  def self.symmetric?
    true
  end

  def self.compare(c1, s1, c2, s2)
  end
end

class CovarianceComparison < CorpusComparison

  # Compute E[(corpus1 - mean(corpus1))(corpus2 - mean(corpus2))]
  def self.compare(corpus1, stats1, corpus2, stats2)
    # Combine key list
    keys = (corpus1.keys + corpus2.keys).uniq

    # Counts
    count = 0
    cov = 0.0

    keys.each{|k|
      cov += (corpus1[k].to_f - stats1[:mean]) * (corpus2[k].to_f - stats2[:mean])
      count += 1
    }

    return cov.to_f / count
  end
end


class PearsonRComparison < CovarianceComparison
  def self.compare(corpus1, stats1, corpus2, stats2)
    cov = super(corpus1, stats1, corpus2, stats2)
    return cov / (stats1[:sd] * stats2[:sd])
  end
end

class KendallTauAComparison < CorpusComparison
  # TODO: rewrite algorithm as nlogn version based
  # on merge sort.
  def self.compare(corpus1, stats1, corpus2, stats2)
    concordant_pairs = 0
    discordant_pairs = 0

    # Combine key list
    keys = (corpus1.keys + corpus2.keys).uniq
    count = keys.length

    keys.each_index{|i|
      ((i+1)..keys.length).each{|j|
        if corpus1[keys[i]] and corpus2[keys[i]] and corpus1[keys[j]] and corpus2[keys[j]] then
          v1i = corpus1[keys[i]]
          v1j = corpus1[keys[j]]
          v2i = corpus2[keys[i]]
          v2j = corpus2[keys[j]]

          if (v1i > v1j) == (v2i > v2j) 
            concordant_pairs += 1
          else
            discordant_pairs += 1
          end
        end
      }
    }

    return (concordant_pairs - discordant_pairs).to_f / (0.5 * count * (count - 1))
  end
end


class KendallTauBComparison < CorpusComparison
  # TODO: rewrite algorithm as nlogn version based
  # on merge sort.
  def self.compare(corpus1, stats1, corpus2, stats2)
    concordant_pairs = 0
    discordant_pairs = 0

    # ties
    t1 = 0
    t2 = 0

    # Combine key list
    keys = (corpus1.keys + corpus2.keys).uniq
    count = keys.length

    keys.each_index{|i|
      ((i+1)..keys.length).each{|j|
        v1i = corpus1[keys[i]]
        v1j = corpus1[keys[j]]
        v2i = corpus2[keys[i]]
        v2j = corpus2[keys[j]]

        if v1i == v1j 
          t1 += 1
        elsif v2i == v2j
          t2 += 1
        elsif (v1i.to_f > v1j.to_f) == (v2i.to_f > v2j.to_f) 
          concordant_pairs += 1
        else
          discordant_pairs += 1
        end
      }
    }

    denom = ( (concordant_pairs + discordant_pairs + t1) * (concordant_pairs + discordant_pairs + t2) ) ** 0.5

    return (concordant_pairs - discordant_pairs).to_f / denom
  end
end


class StuartTauCComparison < CorpusComparison
  # TODO: rewrite algorithm as nlogn version based
  # on merge sort.
  # http://www.unesco.org/webworld/idams/advguide/Chapt4_2.htm
  # https://stat.ethz.ch/pipermail/r-help/2006-September/112806.html
  def self.compare(corpus1, stats1, corpus2, stats2)
    concordant_pairs = 0
    discordant_pairs = 0

    # Combine key list
    keys = (corpus1.keys + corpus2.keys).uniq
    count = keys.length

    keys.each_index{|i|
      keys.each_index{|j|
        if i > j then
          v1i = corpus1[keys[i]]
          v1j = corpus1[keys[j]]
          v2i = corpus2[keys[i]]
          v2j = corpus2[keys[j]]

          # Compare {(v1i, v2i), (v1j, v2j)}
          # FIXME: Check for ties
          if (v1i.to_f > v1j.to_f) == (v2i.to_f > v2j.to_f) 
            concordant_pairs += 1
          else
            discordant_pairs += 1
          end
        end
      }
    }


    m = [corpus1.length, corpus2.length].min.to_f
    denom = ( (count**2 * (m-1)) )

    return (concordant_pairs - discordant_pairs) * ( (2*m) / denom )
  end
end

class GoodmanKruskalGammaComparison < CorpusComparison
  def self.compare(corpus1, stats1, corpus2, stats2)
    concordant_pairs = 0
    discordant_pairs = 0

    # Combine key list
    keys = (corpus1.keys + corpus2.keys).uniq

    keys.each_index{|i|
      keys.each_index{|j|
        if i > j then
          v1i = corpus1[keys[i]]
          v1j = corpus1[keys[j]]
          v2i = corpus2[keys[i]]
          v2j = corpus2[keys[j]]

          # Compare {(v1i, v2i), (v1j, v2j)}
          # Check for ties
          if (v1i.to_f > v1j.to_f) == (v2i.to_f > v2j.to_f) 
            concordant_pairs += 1
          else
            discordant_pairs += 1
          end
        end
      }
    }

    return (concordant_pairs - discordant_pairs).to_f / (concordant_pairs + discordant_pairs)
  end
end


class AsymSomersDComparison < CorpusComparison
  def self.symmetric?
    false
  end

  def self.compare(corpus1, stats1, corpus2, stats2)
    concordant_pairs = 0
    discordant_pairs = 0

    # ties only on corpus 2
    ties = 0

    # Combine key list
    keys = (corpus1.keys + corpus2.keys).uniq

    keys.each_index{|i|
      keys.each_index{|j|
        if i > j then
          v1i = corpus1[keys[i]]
          v1j = corpus1[keys[j]]
          v2i = corpus2[keys[i]]
          v2j = corpus2[keys[j]]

          # Compare {(v1i, v2i), (v1j, v2j)}
          # Check for ties
          if v2i == v2j then
            ties += 1
          elsif (v1i.to_f > v1j.to_f) == (v2i.to_f > v2j.to_f) 
            concordant_pairs += 1
          else
            discordant_pairs += 1
          end
        end
      }
    }  

    return (concordant_pairs + discordant_pairs).to_f / (concordant_pairs + discordant_pairs + ties)
  end
end

class SymSomersDComparison < CorpusComparison
  def self.compare(corpus1, stats1, corpus2, stats2)
    concordant_pairs = 0
    discordant_pairs = 0

    # ties only on corpus 2
    ties1 = 0
    ties2 = 0

    # Combine key list
    keys = (corpus1.keys + corpus2.keys).uniq

    keys.each_index{|i|
      keys.each_index{|j|
        if i > j then
          v1i = corpus1[keys[i]]
          v1j = corpus1[keys[j]]
          v2i = corpus2[keys[i]]
          v2j = corpus2[keys[j]]

          # Compare {(v1i, v2i), (v1j, v2j)}
          # Check for ties
          if v1i == v1j then
            ties1 += 1
          elsif v2i == v2j
            ties2 += 1
          elsif (v1i.to_f > v1j.to_f) == (v2i.to_f > v2j.to_f) 
            concordant_pairs += 1
          else
            discordant_pairs += 1
          end
        end
      }
    } 



    return (concordant_pairs + discordant_pairs).to_f / (0.5 * (concordant_pairs + discordant_pairs + ties1 + concordant_pairs + discordant_pairs + ties2))
  end
end


class ChiSquareComparison < CorpusComparison
  def self.symmetric?
    false
  end

  # NB: corpus 1 is the "expected value" corpus, null hyp, etc.
  def self.compare(corpus1, stats1, corpus2, stats2)
    # Combine key list
    keys = (corpus1.keys + corpus2.keys).uniq

    chisq = 0.0

    keys.each{|k|
      expected = corpus1[k].to_f + 0.001  # TODO: check this adjustment
      freq = corpus2[k].to_f


      chisq += ((freq - expected)**2) / expected
    }

    puts "chisq: #{chisq}"# #{((freq - expected)**2) / expected}"
    # TODO: normalise and look up using chisq distribution
    # tables.
    return chisq
  end
end

class LikelihoodRatioComparison < CorpusComparison
  def self.symmetric?
    false
  end

  # NB: corpus 1 is the "expected value" corpus, null hyp, etc.
  def self.compare(corpus1, stats1, corpus2, stats2)
    # Combine key list
    keys = (corpus1.keys + corpus2.keys).uniq

    lr = 0

    keys.each{|k|
      expected = corpus1[k].to_f + 0.001  # TODO: check this adjustment
      freq = corpus2[k].to_f

      lr += freq * Math.log( freq / expected ) if freq > 0
    }

    return lr*2
  end
end

class MantelHaenszelComparison < CorpusComparison
  # NB: corpus 1 is the "expected value" corpus, null hyp, etc.
  # 
  # FIXME: this needs checking to ensure the value of W makes sense.
  def self.compare(corpus1, stats1, corpus2, stats2)
    r = PearsonRComparison.compare(corpus1, stats1, corpus2, stats2)
    return ((stats1[:sum] + stats2[:sum]) - 1) * (r**2)
  end
end


class SpearmanComparison < CorpusComparison
  def self.compare(corpus1, stats1, corpus2, stats2)
    r1 = rank_corpus(corpus1)
    s1 = stats(r1, {})
    r2 = rank_corpus(corpus2)
    s2 = stats(r2, {})

    return PearsonRComparison.compare(r1, s1, r2, s2)
  end

  private
  def self.rank_corpus(corpus)
    ranked = {}
    c = -1
    corpus.sort_by{|k,v| v}.each{|k, v|
      ranked[k] = (c += 1)
    }
    return ranked
  end
end

# square only
class CohenKappaComparison < CorpusComparison
  def self.compare(corpus1, stats1, corpus2, stats2)
  end
end


class TschuprowTComparison < CorpusComparison
  def self.compare(corpus1, stats1, corpus2, stats2)
  end
end

# FIXME: this is just wrong, should be 1 when variables are identical
class CramerVComparison < CorpusComparison
  def self.compare(corpus1, stats1, corpus2, stats2)
    chisq = ChiSquareComparison.compare(corpus1, stats1, corpus2, stats2)

    k = [corpus1, corpus2].map{|x| x.length}.min
    n = stats1[:sum] + stats2[:sum]
      # corpus1.length * corpus2.length # or this: (corpus1.keys + corpus2.keys).uniq ^ 2

    puts "cramerv: #{chisq}, #{k}, #{n}"

    return ( chisq / ( n * (k - 1) ) ) ** 0.5
  end
end


class GoodmanKruskalLambdaComparison < CorpusComparison
  def self.compare(corpus1, stats1, corpus2, stats2)
  end
end


class KilgarriffComparison < CorpusComparison
  def self.compare(corpus1, stats1, corpus2, stats2)
  end
end

# print out the full set for two corpora, along with correlation (cov(i,j)/sdi*sdj)
def similarity_matrix(comparison, corpora, stats)

  cov =  []
  corpora.length.times{|i|  cov[i] = [] }

  corpora.length.times{|i|
    corpora.length.times{|j|
      if not comparison.symmetric? or (i >= j) then
        cov[i][j] = comparison.compare(corpora[i], stats[i], corpora[j], stats[j])


        # cov[i][j] = covariance(corpora[i], stats[i], corpora[j], stats[j])
        # cor[i][j] = cov[i][j] / (stats[i][:sd] * stats[j][:sd])

        # puts "COV(#{i},#{j}) = #{covariance(corpora[i], stats[i], corpora[j], stats[j])}"
        # puts "PMCC(#{i},#{j}) = #{(covariance(corpora[i], stats[i], corpora[j], stats[j])) / (stats[i][:sd] * stats[j][:sd])}"
      end
    }
  }

  return cov
end

def print_matrix(data)
  puts "\t#{(1..data.length).to_a.join("\t")}"

  data.each_index{|r|
    print "#{r+1}\t" 

    data[r].each_index{|c|
      d = data[r][c].to_f
      print "#{(d<0)?'':' '}#{d.round(2)} \t"
    }
    print "\n"
  }
end

def test(klass, corpora, stats)
  puts "\n"
  puts "#{klass}"
  puts "-" * klass.to_s.length
  print_matrix(similarity_matrix(klass, corpora, stats))
end


# Compute mean, variance
stats_one = stats(corpus_one, stats_one)
stats_two = stats(corpus_two, stats_two)

puts "1: #{stats_one}"
puts "2: #{stats_two}"

# test(CovarianceComparison, [corpus_one, corpus_two], [stats_one, stats_two])
test(PearsonRComparison, [corpus_one, corpus_two], [stats_one, stats_two])
# test(KendallTauAComparison, [corpus_one, corpus_two], [stats_one, stats_two])
# test(KendallTauBComparison, [corpus_one, corpus_two], [stats_one, stats_two])
# test(StuartTauCComparison, [corpus_one, corpus_two], [stats_one, stats_two])
# test(GoodmanKruskalGammaComparison, [corpus_one, corpus_two], [stats_one, stats_two])
# test(AsymSomersDComparison, [corpus_one, corpus_two], [stats_one, stats_two])
# test(SymSomersDComparison, [corpus_one, corpus_two], [stats_one, stats_two])
test(SpearmanComparison, [corpus_one, corpus_two], [stats_one, stats_two])
# test(CohenKappaComparison, [corpus_one, corpus_two], [stats_one, stats_two])
test(ChiSquareComparison, [corpus_one, corpus_two], [stats_one, stats_two])
# test(LikelihoodRatioComparison, [corpus_one, corpus_two], [stats_one, stats_two])
# test(MantelHaenszelComparison, [corpus_one, corpus_two], [stats_one, stats_two])
# test(CramerVComparison, [corpus_one, corpus_two], [stats_one, stats_two])



# --
# test(Eta, [corpus_one, corpus_two], [stats_one, stats_two])
# test(Eta, [corpus_one, corpus_two], [stats_one, stats_two])
# test(Eta, [corpus_one, corpus_two], [stats_one, stats_two])
# test(Eta, [corpus_one, corpus_two], [stats_one, stats_two])
# test(CorrelationComparison, [corpus_one, corpus_two], [stats_one, stats_two])
# test(KilgarriffComparison, [corpus_one, corpus_two], [stats_one, stats_two])
#
#
#
# from http://salises.mona.uwi.edu/sa63c/Crosstabs%20Measures%20for%20Nominal%20Data.htm :
# The statistics options for nominal measures include: the contingency coefficient, phi, Cramér's V, lambda, and the uncertainty coefficient. Three of these measures, phi, Cramer's V, and the contingency coefficient, are based on the chi-square itself
# Statistics for ordinal measures include gamma, Somers' d, Kendall's tau-b, and Kendall's tau-c.
# There is one statistic, eta, designed for the case when one one of the measures is nominal and the other is interval.
# The correlation option yields both Spearman's rank order correlation (used for ordinal data) and Pearson's product moment correlation (used for interval data).
# Other statistics include kappa, used as a measure of inter-rater reliability, a risk analysis test, and McNemar's test, used as test of change in repeated measures designs.
# The statistics discussed in this set of notes include chi-square and the nominal statistics (contingency coefficient, phi, Cramer's V, lambda, and the uncertainty coefficient).
#
