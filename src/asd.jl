using Dates

## --------------------------------struct Definition------------------------------------------

mutable struct Header
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
    piezoExtensionX   ::Float64    # Xピエゾ伸び係数[nm/V]
    piezoExtensionY   ::Float64    # Yピエゾ伸び係数[nm/V]
    piezoExtensionZ   ::Float64    # Zピエゾ伸び係数[nm/V]
    piezoGainZ        ::Float64
    adRange           ::Float64    # AD電圧レンジ
    AdResolution      ::Float64    # AD分解能
    isAveraged        ::Bool       # 移動平均化フラグ（tureで移動平均）
    averageWindow     ::Int64      # 1ピクセルに使える最大データ数(最低１)
    day               ::DateTime
    roundingDegreeX   ::Int64
    roundingDegreeY   ::Int64
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
    textEncodingCodepage ::Int64

    Header() = new()
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
    header   ::Header
    frames   ::Array{Frame, 1}
end

## --------------------------------struct Definition------------------------------------------
## -----------------------------------read helper---------------------------------------------

function getDataType(data)
    if     data == 0x5054 return "topography"
    elseif data == 0x5245 return "error"
    elseif data == 0x4850 return "phase"
    else return "none" end
end

function getAdRange(data)
    # ADレンジ定数（ユニポーラ0～1V）(使われていないらしい)
    if     data == 0x00000001 @assert false "invalid AdRange"
    # ADレンジ定数（ユニポーラ0～2.5V）(使われていないらしい)
    elseif data == 0x00000002 @assert false "invalid AdRange"
    # ADレンジ定数（ユニポーラ0～5V）(使われていないらしい)
    elseif data == 0x00000004 @assert false "invalid AdRange"
    # ADレンジ定数（バイポーラ±1V）
    elseif data == 0x00010000 return 2.0
    # ADレンジ定数（バイポーラ±2.5V）
    elseif data == 0x00020000 return 5.0
    # ADレンジ定数（バイポーラ±5V）
    elseif data == 0x00040000 return 10.0
    # ADレンジ定数（バイポーラ±80V, データを編集した場合に仮想的にこれを使う。実際にバイポーラ80VでAD変換したわけではない。また、分解能は16ビットにする）
    elseif data == 0x00080000 return 160.0
    # 何も当てはまらない
        else @assert false "invalid AdRange" end
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

## -----------------------------------read helper---------------------------------------------

function readDateV0(io::IOStream)
    year              = Int64(read(io, Int16))
    month             = Int64(read(io, UInt8))
    day               = Int64(read(io, UInt8))
    hour              = Int64(read(io, UInt8))
    minute            = Int64(read(io, UInt8))
    second            = Int64(read(io, UInt8))
    return DateTime(year, month, day, hour, minute, second)
end

function readDateV1(io::IOStream)
    year              = Int64(read(io, Int32))
    month             = Int64(read(io, Int32))
    day               = Int64(read(io, Int32))
    hour              = Int64(read(io, Int32))
    minute            = Int64(read(io, Int32))
    second            = Int64(read(io, Int32))
    return DateTime(year, month, day, hour, minute, second)
end

function readHeaderV0(io::IOStream)
    header = Header()
    header.fileVersion       = 0
    header.dataType1ch       = getDataType(read(io, Int16))
    header.dataType2ch       = getDataType(read(io, Int16))
    header.fileHeaderSize    = Int64(read(io, Int32))
    header.frameHeaderSize   = Int64(read(io, Int32))
    header.operatorNameSize  = Int64(read(io, Int32))
    header.commentOffsetSize = Int64(read(io, Int32))
    header.commentSize       = Int64(read(io, Int32))
    header.pixelX            = Int64(read(io, Int16))
    header.pixelY            = Int64(read(io, Int16))
    header.scanningRangeX    = Int64(read(io, Int16))
    header.scanningRangeY    = Int64(read(io, Int16))
    header.frameRate         = Float64(read(io, Float32))
    header.piezoExtensionZ   = Float64(read(io, Float32))
    header.piezoGainZ        = Float64(read(io, Float32))
    header.adRange           = getAdRange(read(io, UInt32))
    header.AdResolution      = (2.0)^Int64(read(io, Int32))
    header.isAveraged        = read(io, Bool)
    header.averageWindow     = Int64(read(io, Int16))
    legacy                   = Int64(read(io, Int32)) # ダミー
    header.day               = readDateV0(io)
    header.roundingDegreeX   = Int64(read(io, UInt8))
    header.maxRangeX         = Float64(read(io, Float32))
    header.maxRangeY         = Float64(read(io, Float32))
    header.booked1           = Int64(read(io, Int32))
    header.booked2           = Int64(read(io, Int32))
    header.booked3           = Int64(read(io, Int32))
    header.initFrame         = Int64(read(io, Int32))
    header.numFrames         = Int64(read(io, Int32))
    header.machineId         = Int64(read(io, Int32))
    header.fileId            = Int64(read(io, Int16))
    header.operatorName      = getOperatorName(io, header.operatorNameSize)
    header.sensorSensitivity = Float64(read(io, Float32))
    header.phaseSensitivity  = Float64(read(io, Float32))
    header.scannigDirection  = Int64(read(io, Int32))
    header.comment           = getComment(io, header.commentOffsetSize, header.commentSize)

    return header
end

function readHeaderV1(io::IOStream)
    header = Header()
    header.fileVersion       = 1
    header.fileHeaderSize    = Int64(read(io, Int32))
    header.frameHeaderSize   = Int64(read(io, Int32))
    header.textEncodingCodepage = Int64(read(io, Int32))
    header.operatorNameSize  = Int64(read(io, Int32))
    header.commentSize       = Int64(read(io, Int32))
    header.dataType1ch       = getDataType(read(io, Int32))
    header.dataType2ch       = getDataType(read(io, Int32))
    header.initFrame         = Int64(read(io, Int32))
    header.numFrames         = Int64(read(io, Int32))
    header.scannigDirection  = Int64(read(io, Int32))
    header.fileId            = Int64(read(io, Int32))
    header.pixelX            = Int64(read(io, Int32))
    header.pixelY            = Int64(read(io, Int32))
    header.scanningRangeX    = Int64(read(io, Int32))
    header.scanningRangeY    = Int64(read(io, Int32))
    header.isAveraged        = read(io, Bool)
    header.averageWindow     = Int64(read(io, Int32))
    header.day               = readDateV1(io)
    header.roundingDegreeX   = Int64(read(io, Int32))
    header.roundingDegreeY   = Int64(read(io, Int32))
    header.frameRate         = Float64(read(io, Float32))
    header.sensorSensitivity = Float64(read(io, Float32))
    header.phaseSensitivity  = Float64(read(io, Float32))
    header.commentOffsetSize = Int64(read(io, Int32))
    header.booked1           = Int64(read(io, Int32))
    header.booked2           = Int64(read(io, Int32))
    header.booked3           = Int64(read(io, Int32))
    header.machineId         = Int64(read(io, Int32))
    header.adRange           = getAdRange(read(io, Int32))
    header.AdResolution      = (2.0)^Int64(read(io, Int32))
    header.maxRangeX         = Float64(read(io, Float32))
    header.maxRangeY         = Float64(read(io, Float32))
    header.piezoExtensionX   = Float64(read(io, Float32))
    header.piezoExtensionY   = Float64(read(io, Float32))
    header.piezoExtensionZ   = Float64(read(io, Float32))
    header.piezoGainZ        = Float64(read(io, Float32))
    header.operatorName      = getOperatorName(io, header.operatorNameSize)
    header.comment           = getComment(io, 0, header.commentSize)

    return header
end

"""
ADのバイナリ-物理量の変換公式は
物理量 = ボードのAD変換レンジ*(サンプリングバイナリデータ)/2^(分解能) - ボードのAD変換レンジの半分。
また、PIDの信号を高さ情報として取り込んでいるので、

バイナリデータが大きい = PIDの出力電圧が高い = 試料に対して押し込んでいる = 高さが低い

という関係になる。
したがって、高さの最大最小は

高さ最値小 = バイナリ最大値から計算した高さ
高さ最大値 = バイナリ最小値から計算した高さ
"""
function binaryToPhysicalQuantity!(data, header, chanelType, unit)
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
        @assert false "invalid chanelType"
    end

    if unit == "angstrom"
        # nm -> angstrom(要議論)
        unitConversion = 10.0
    else
        # nm
        unitConversion = 1.0
    end

    for y in 1:header.pixelY, x in 1:header.pixelX
        data[y, x] = (adUiniRange - data[y, x] * cc) * multiplier * unitConversion
    end
end

function readFrameHeader(io::IOStream, header::Header)
    io_pos = position(io)
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
    # seek(io, io_pos + header.frameHeaderSize)

    return FrameHeader(number, maxData, minData, offsetX, offsetY, tiltX, tiltY, isStimulated, booked1, booked2, booked3, booked4)
end

function readImage(io::IOStream, header, chanelType, unit)
    data = zeros(header.pixelY, header.pixelX)
    # TODO: 平均化回数が1じゃないことがあるらしい(今は不必要？)
    for y in 1:header.pixelY, x in 1:header.pixelX
        data[y, x] = Int64(read(io, Int16))
    end
    binaryToPhysicalQuantity!(data, header, chanelType, unit)
    return data
end

function readFrameRangeChack(readFrameRange, header)
    if readFrameRange == nothing
        readFrameRange = (1, header.numFrames)
    end
    if readFrameRange[1] < 1 || readFrameRange[2] > header.numFrames
        @printf "readFrameRange is out of Frame size"
    end
    if readFrameRange[1] > readFrameRange[2]
        @printf "readFrameRange is invalid"
    end
    readFrameRange = (max(min(readFrameRange[1], header.numFrames), 1), min(max(readFrameRange[2], 1), header.numFrames))
    return readFrameRange
end

function makeExpandedAsd(header, readFrameRange, frameHeaders, datas, subDatas, maxFrameSize)
    frames = []
    offsetY_min = Int64(2^63 - 1)
    offsetY_max = Int64(-2^63)
    offsetX_min = Int64(2^63 - 1)
    offsetX_max = Int64(-2^63)
    for i in 1:(readFrameRange[2] - readFrameRange[1] + 1)
        offsetY_min = min(offsetY_min, frameHeaders[i].offsetY)
        offsetY_max = max(offsetY_max, frameHeaders[i].offsetY)
        offsetX_min = min(offsetX_min, frameHeaders[i].offsetX)
        offsetX_max = max(offsetX_max, frameHeaders[i].offsetX)
    end
            
    expandedY = offsetY_max - offsetY_min
    expandedX = offsetX_max - offsetX_min
    
    if ((header.pixelY + expandedY) * (header.pixelX + expandedX)) > maxFrameSize
        @printf "ERROR: frame size will be %d, it's too big" (header.pixelY + expandedY) * (header.pixelX + expandedX)
        return Asd(header, [])
    end
            
    for i in 1:(readFrameRange[2] - readFrameRange[1] + 1)
        y_min = frameHeaders[i].offsetY - offsetY_min + 1
        y_max = y_min + header.pixelY - 1
        x_min = frameHeaders[i].offsetX - offsetX_min + 1
        x_max = x_min + header.pixelX - 1
        data = zeros(header.pixelY + expandedY, header.pixelX + expandedX)
        data[y_min:y_max, x_min:x_max] = datas[i]
        subData = zeros(header.pixelY + expandedY, header.pixelX + expandedX)
        subData[y_min:y_max, x_min:x_max] = subDatas[i]
        push!(frames, Frame(frameHeaders[i], data, subData))
    end
            
    return Asd(Header(), frames)
end

"""
    readasd(filePath; readFrameRange=nothing, translationSetting=nothing, unit="angstrom", maxFrameSize=1000000) -> output::Asd

Read a specified ASD file. The user can specify a range of frames to read by `readFrameRange` option.
readFrameRange is a tuple of (startFrame, endFrame). The user can specify a translation setting by `translationSetting` option.
If translationSetting is "expansion", images are expanded so as to cancel out the translational offset (moving of the focus) of each frame.
The unit of output data is specified by `unit` option. If unit is "angstrom" (default), the unit is converted to angstrom.
Otherwise, the unit is kept as the original unit (nanometer by default).

Return a `Asd` object as output. The Asd objects contains header and (multiple) frame(s) information.

# Example
```julia-repl
julia> asd = mdload("ak.asd")
julia> asd.header
julia> asd.frames[1]
julia> using Plots
julia> heatmap(asd.frames[1].data)
```
"""
function readasd(filePath; readFrameRange = nothing, translationSetting = nothing, unit = "angstrom", maxFrameSize = 1000000)
    open(filePath, "r") do io
        # io_pos = position(io)
        headerVersion = Int64(read(io, Int32))
        header = Header()
        if headerVersion == 0
            header = readHeaderV0(io)
        elseif headerVersion == 1
            header = readHeaderV1(io)
        else
            @assert false "can't read v2 file"
        end
        
        header.fileHeaderSize = position(io)
        readFrameHeader(io, header)
        readImage(io, header, header.dataType1ch, unit)
        oneFrameSize = position(io) - header.fileHeaderSize
        seek(io, header.fileHeaderSize)

        frameHeaders = []
        subFrameHeaders = []
        datas = []
        subDatas = []
        
        readFrameRange = readFrameRangeChack(readFrameRange, header)
        
        for i in readFrameRange[1]:readFrameRange[2]
            seek(io, header.fileHeaderSize + oneFrameSize * (i - 1))
            push!(frameHeaders, readFrameHeader(io, header))
            push!(datas, readImage(io, header, header.dataType1ch, unit))
        end

        for i in readFrameRange[1]:readFrameRange[2]
            if header.dataType2ch == "none"
                push!(subDatas, datas[i - readFrameRange[1] + 1])
                continue
            end
            seek(io, header.fileHeaderSize + oneFrameSize * (header.numFrames + i - 1))
            push!(subFrameHeaders, readFrameHeader(io, header))
            push!(subDatas, readImage(io, header, header.dataType1ch, unit))
        end

        if translationSetting == "expansion"
            return makeExpandedAsd(header, readFrameRange, frameHeaders, datas, subDatas, maxFrameSize)
        end
        
        frames = []
        for i in 1:(readFrameRange[2] - readFrameRange[1] + 1)
            push!(frames, Frame(frameHeaders[i], datas[i], subDatas[i]))
        end

        return Asd(header, frames)
    end
end
