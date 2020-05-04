using Dates

abstract type AbstractHeader end

struct HeaderV0 <: AbstractHeader
    fileVersion       ::Int64
    dataType1ch       ::String     # 1chのデータ種別
    dataType2ch       ::String     # 2chのデータ種別
    fileHeaderSize    ::Int64      # ファイルヘッダサイズ（コメントのサイズ含む。）
    frameHeaderSize   ::Int64      # フレームヘッダサイズ（コメントのサイズ含む。）
    operatorNameSize  ::Int64
    commentOffsetSize ::Int64
    commentSize       ::Int64
    pixelX            ::Int64
    pixelY            ::Int64
    scanningRangeX    ::Int64
    scanningRangeY    ::Int64
    frameRate         ::Float64
    piezoExtensionZ   ::Float64    # Zピエゾ伸び係数[nm/V]
    piezoGainZ        ::Float64
    adRange           ::Float64    # AD電圧レンジ
    AdResolution      ::Float64    # AD分解能
    isAveraged        ::Bool       # 移動平均化フラグ（tureで移動平均）
    averageWindow     ::Int64      # 1ピクセルに使える最大データ数(最低１)
    day               ::DateTime
    roundingDegree    ::Int64
    maxRangeX         ::Float64
    maxRangeY         ::Float64
    booked1           ::Int64
    booked2           ::Int64
    booked3           ::Int64
    initFrame         ::Int64
    numFrames         ::Int64
    machineId         ::Int64
    fileId            ::Int64
    operatorName      ::Array{Int, 1}
    sensorSensitivity ::Float64       # センサー感度
    phaseSensitivity  ::Float64       # 位相感度
    scannigDirection  ::Int64         # データ取得ステータスコード
    comment           ::Array{Int, 1}
end

struct FrameHeader
    number       ::Int64
    maxData      ::Int64
    minData      ::Int64
    offsetX      ::Int64
    offsetY      ::Int64
    tiltX        ::Float64
    tiltY        ::Float64
    isStimulated ::Bool
    booked1      ::Int64
    booked2      ::Int64
    booked3      ::Int64
    booked4      ::Int64
end

struct Frame
    header   ::FrameHeader
    data     ::Array{Float64, 2}
    subData  ::Union{Array{Float64, 2}, Nothing}
end

# 現在単位はÅ
struct Asd
    header   ::AbstractHeader
    frames   ::Array{Frame, 1}
end

function getDataType(data)
    if     data == 0x5054 return "topography"
    elseif data == 0x5245 return "error"
    elseif data == 0x4850 return "phase"
    else return "none" end
end

function getAdRange(data)
    # ADレンジ定数（ユニポーラ0～1V）(使われていないらしい)
    if     data == 0x00000001 @assert false
    # ADレンジ定数（ユニポーラ0～2.5V）(使われていないらしい)
    elseif data == 0x00000002 @assert false
    # ADレンジ定数（ユニポーラ0～5V）(使われていないらしい)
    elseif data == 0x00000004 @assert false
    # ADレンジ定数（バイポーラ±1V）
    elseif data == 0x00010000 return 2.0
    # ADレンジ定数（バイポーラ±2.5V）
    elseif data == 0x00020000 return 5.0
    # ADレンジ定数（バイポーラ±5V）
    elseif data == 0x00040000 return 10.0
    # ADレンジ定数（バイポーラ±80V, データを編集した場合に仮想的にこれを使う。実際にバイポーラ80VでAD変換したわけではない。また、分解能は16ビットにする）
    elseif data == 0x00080000 return 160.0
    # 何も当てはまらない
    else @assert false end
    return nothing
end

function getOperatorName(io, operatorNameSize)
    ret = []
    for i in 1:operatorNameSize
        push!(ret, Int64(read(io, UInt8)))
    end
    return ret
end

function getComment(io, commentOffsetSize, commentSize)
    for i in 1:commentOffsetSize
        read(io, Bool)
    end
    ret = []
    for i in 1:commentSize
        push!(ret, Int64(read(io, UInt8)))
    end
    return ret
end

function readDate(io::IOStream)
    year              = Int64(read(io, Int16))
    month             = Int64(read(io, UInt8))
    day               = Int64(read(io, UInt8))
    hour              = Int64(read(io, UInt8))
    minute            = Int64(read(io, UInt8))
    second            = Int64(read(io, UInt8))
    return DateTime(year, month, day, hour, minute, second)
end

function readHeaderV0(io::IOStream)
    fileVersion       = 0
    dataType1ch       = getDataType(read(io, Int16))
    dataType2ch       = getDataType(read(io, Int16))
    fileHeaderSize    = Int64(read(io, Int32))
    frameHeaderSize   = Int64(read(io, Int32))
    operatorNameSize  = Int64(read(io, Int32))
    commentOffsetSize = Int64(read(io, Int32))
    commentSize       = Int64(read(io, Int32))
    pixelX            = Int64(read(io, Int16))
    pixelY            = Int64(read(io, Int16))
    scanningRangeX    = Int64(read(io, Int16))
    scanningRangeY    = Int64(read(io, Int16))
    frameRate         = Float64(read(io, Float32))
    piezoExtensionZ   = Float64(read(io, Float32))
    piezoGainZ        = Float64(read(io, Float32))
    adRange           = getAdRange(read(io, UInt32))
    AdResolution      = (2.0)^Int64(read(io, Int32))
    isAveraged        = read(io, Bool)
    averageWindow     = Int64(read(io, Int32))
    legacy            = Int64(read(io, Int16)) # ダミー
    day               = readDate(io)
    roundingDegree    = Int64(read(io, UInt8))
    maxRangeX         = Float64(read(io, Float32))
    maxRangeY         = Float64(read(io, Float32))
    booked1           = Int64(read(io, Int32))
    booked2           = Int64(read(io, Int32))
    booked3           = Int64(read(io, Int32))
    initFrame         = Int64(read(io, Int32))
    numFrames         = Int64(read(io, Int32))
    machineId         = Int64(read(io, Int32))
    fileId            = Int64(read(io, Int16))
    operatorName      = getOperatorName(io, operatorNameSize)
    sensorSensitivity = Float64(read(io, Float32))
    phaseSensitivity  = Float64(read(io, Float32))
    scannigDirection  = Int64(read(io, Int32))
    comment           = getComment(io, commentOffsetSize, commentSize)

    return HeaderV0(fileVersion, dataType1ch, dataType2ch, fileHeaderSize, frameHeaderSize, operatorNameSize, commentOffsetSize, commentSize, pixelX, pixelY, scanningRangeX, scanningRangeY, frameRate, piezoExtensionZ, piezoGainZ, adRange, AdResolution, isAveraged, averageWindow, day, roundingDegree, maxRangeX, maxRangeY, booked1, booked2, booked3, initFrame, numFrames, machineId, fileId, operatorName, sensorSensitivity, phaseSensitivity, scannigDirection, comment)
end

"""
ADのバイナリ－物理量の変換公式は
物理量 = ボードのAD変換レンジ*(サンプリングバイナリデータ)/2^(分解能) - ボードのAD変換レンジの半分。
また、PIDの信号を高さ情報として取り込んでいるので、

バイナリデータが大きい　＝　PIDの出力電圧が高い　＝　試料に対して押し込んでいる　＝　高さが低い

という関係になる。
したがって、高さの最大最小は

高さ最値小 = バイナリ最大値から計算した高さ
高さ最大値 = バイナリ最小値から計算した高さ
"""
function binaryToPhysicalQuantity(data, header, chanelType)
    # 電圧データに変換
    cc = header.adRange / header.AdResolution
    adUiniRange = header.adRange / 2

    # 電圧データを高さor位相データに変換するための乗数。ErrorとPhaseの場合はTopo像とは逆符号になる点に注意。
    multiplier = Float64(0)
    if chanelType == "topography"
        multiplier = header.piezoGainZ * header.piezoExtensionZ
    elseif chanelType == "error"
        multiplier = -1.0 * header.sensorSensitivity
    elseif chanelType == "phase"
        multiplier = if header.phaseSensitivity != 0 phaseSensitivity else -1.0 end
    else
        # ここには来ない
        @assert false
    end

    # nm -> angstrom(要議論)
    unitConversion = 10

    for y in 1:header.pixelY, x in 1:header.pixelX
        data[y, x] = (adUiniRange - data[y, x] * cc) * multiplier * unitConversion
    end
end

function readFrameHeader(io::IOStream)
    number       = Int64(read(io, Int32))
    maxData      = Int64(read(io, Int16))
    minData      = Int64(read(io, Int16))
    offsetX      = Int64(read(io, Int16))
    offsetY      = Int64(read(io, Int16))
    tiltX        = Float64(read(io, Float32))
    tiltY        = Float64(read(io, Float32))
    isStimulated = read(io, Bool)
    booked1      = Int64(read(io, Int8))
    booked2      = Int64(read(io, Int16))
    booked3      = Int64(read(io, Int32))
    booked4      = Int64(read(io, Int32))

    return FrameHeader(number, maxData, minData, offsetX, offsetY, tiltX, tiltY, isStimulated, booked1, booked2, booked3, booked4)
end

function readImage(io::IOStream, header, chanelType)
    data = zeros(header.pixelY, header.pixelX)
    # TODO: 平均化回数が1じゃないことがあるらしい(今は不必要？)
    for y in 1:header.pixelY, x in 1:header.pixelX
        data[y, x] = Int64(read(io, Int16))
    end
    binaryToPhysicalQuantity(data, header, chanelType)
    return data
end

function readFrame(io::IOStream, header)
    frameHeader = readFrameHeader(io)

    data = readImage(io, header, header.dataType1ch)

    if header.dataType2ch == "none"
        return Frame(frameHeader, data, nothing)
    end

    subData = readImage(io, header, header.dataType2ch)

    return Frame(frameHeader, data, subData)
end

function readasd(filePath)
    open(filePath, "r") do io
        headerVersion = Int64(read(io, Int32))

        if headerVersion == 0
            header = readHeaderV0(io)
            frames = []
            for i in 1:header.numFrames
                push!(frames, readFrame(io, header))
            end

            return Asd(header, frames)
        else
            # TODO
            @assert false
        end
    end
end
