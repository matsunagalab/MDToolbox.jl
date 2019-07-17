# 多次元対応？

# 確率密度関数を求める
function calc_densities(p, data, sigma)
    density = []
    for i = 1:length(p)
        push!(density, calc_density(p[i], data, sigma))
    end
    
    return density
end

# 各点の密度を求める
function calc_density(p, data, sigma)
    density = 0.0
    for i = 1:length(data)
        
        # 各次元のガウス値を掛け合わせる
        gaussian = 1.0
        for j = 1:length(data[i])
            gaussian *= get_gaussian(p[j], data[i][j], sigma)
        end
        
        density += gaussian
    end
    
    density /= length(data)
    return density
end

# ガウス関数の値を求める
function get_gaussian(dx, mean, sigma)
    return exp((-(mean - dx)^2)/(2*(sigma^2)))/(sqrt(2.0*pi)*sigma)
end
