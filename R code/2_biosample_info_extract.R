
setwd("D:/02_pku/01_research/06_predatory_bacteria_virome/01_phage_genome_info")
library(tidyr)
library(dplyr)
library(stringr)

# 读取文件内容（假设文件名为 "biosample_result.txt"，请根据实际情况调整）
# 如果内容已在 R 字符向量中，可跳过此行
lines <- readLines("biosample_result.txt", warn = FALSE)

# 将多行合并为一个长字符串，并按空行分割为条目块
full_text <- paste(lines, collapse = "\n")
# 用两个以上换行符分割（条目之间有空行）
blocks <- unlist(strsplit(full_text, "\n\n+"))

# 解析经纬度组合字符串（例如 "39.283 N 75.3633 W"）
parse_lat_lon_str <- function(s) {
  if (is.na(s) || s == "" || grepl("not applicable|not collected|missing: third party data", s, ignore.case = TRUE)) {
    return(c(NA_real_, NA_real_))
  }
  # 替换逗号为点（欧洲小数格式）
  s <- gsub(",", ".", s)
  # 提取所有数字（包括小数点）
  nums <- as.numeric(unlist(regmatches(s, gregexpr("\\d+\\.?\\d*", s))))
  # 提取方向字母（N/S/E/W）
  dirs <- unlist(regmatches(s, gregexpr("[NSEW]", s)))
  if (length(nums) < 2) {
    return(c(NA_real_, NA_real_))
  }
  lat <- nums[1]
  lon <- nums[2]
  # 如果方向存在，根据方向调整符号
  if (length(dirs) >= 2) {
    if (toupper(dirs[1]) == "S") lat <- -lat
    if (toupper(dirs[2]) == "W") lon <- -lon
  }
  # 如果没有方向，但数值可能包含负号（例如经度 -95.15），方向已隐含在符号中
  return(c(lat, lon))
}

# 提取单个字段值（如 "/geographic location (latitude)="..."）
extract_attr <- function(block, key) {
  # 转义括号和空格，使用正则匹配 /key="..."（允许空格）
  pattern <- paste0("/", gsub(" ", "\\\\s+", gsub("\\(", "\\\\(", gsub("\\)", "\\\\)", key))), "=\"([^\"]*)\"")
  m <- regmatches(block, regexec(pattern, block))
  if (length(m) > 0 && length(m[[1]]) > 1) {
    return(m[[1]][2])
  } else {
    return(NA_character_)
  }
}

# 处理每个块
result <- lapply(blocks, function(block) {
  # 提取 BioSample
  bio_match <- regmatches(block, regexec("BioSample:\\s*(SAM[NE]A?\\d+)", block))
  biosample <- if (length(bio_match) > 0 && length(bio_match[[1]]) > 1) bio_match[[1]][2] else NA_character_
  
  # 提取组合经纬度字段
  ll_str <- extract_attr(block, "latitude and longitude")
  if (!is.na(ll_str)) {
    coords <- parse_lat_lon_str(ll_str)
    lat <- coords[1]
    lon <- coords[2]
  } else {
    # 尝试分别提取纬度和经度
    lat_str <- extract_attr(block, "geographic location (latitude)")
    lon_str <- extract_attr(block, "geographic location (longitude)")
    # 处理单个值
    parse_single <- function(x) {
      if (is.na(x)) return(NA_real_)
      if (grepl("not applicable|not collected|missing: third party data", x, ignore.case = TRUE)) {
        return(NA_real_)
      }
      x <- gsub(",", ".", x)  # 替换逗号为点
      return(as.numeric(x))
    }
    lat <- parse_single(lat_str)
    lon <- parse_single(lon_str)
  }
  
  data.frame(BioSample = biosample, latitude = lat, longitude = lon, stringsAsFactors = FALSE)
})

# 合并所有结果
df <- do.call(rbind, result)
# 去除可能因空块产生的 NA 行
df <- df[!is.na(df$BioSample), ]

# 查看结果
print(df)

# 可选：保存为 CSV
# write.csv(df, "biosample_coords.csv", row.names = FALSE)
