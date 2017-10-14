import React from 'react'
import PropTypes from 'prop-types'

import style from './style.css'

function attachProps (prevalenceData, wildType) {
  return prevalenceData.map(({position, aminoAcid, ...props}) => ({
    position,
    aminoAcid,
    ...props,
    isWildType: wildType[position - 1] === aminoAcid
  }))
}

function prevalenceDataFilter (
  prevalenceData, requiredSubtype, minimalPercent) {
  return prevalenceData.filter(({subtype, percent, aminoAcid}) => (
    requiredSubtype === subtype &&
    percent > minimalPercent &&
    ['ins', 'del', 'X'].indexOf(aminoAcid) === -1))
}

function groupPrevalenceDataByPosition (prevalenceData) {
  const result = Object.entries(
    prevalenceData.reduce((acc, prev) => {
      const {position} = prev
      acc[position] = acc[position] || []
      acc[position].push(prev)
      return acc
    }, {})
  )
  return result
    .sort(([a], [b]) => a - b)
    .map(([pos, prevalence]) => [
      pos,
      prevalence.sort((a, b) => b.percent - a.percent)
    ])
}

function makeChunks (array, chunkSize) {
  const result = []
  while (array.length > 0) {
    result.push(array.splice(0, chunkSize))
  }
  return result
}

function smartRound (number) {
  if (number > 1) {
    return Math.round(number)
  } else {
    return Math.round(number * 10) / 10
  }
}

function makeIndelsMap (prevalenceData, requiredSubtype, minimalPercent) {
  return prevalenceData
    .filter(({subtype, percent, aminoAcid}) => (
      requiredSubtype === subtype &&
      percent > minimalPercent &&
      ['ins', 'del'].indexOf(aminoAcid) > -1))
    .reduce((acc, {position}) => {
      acc[position] = true
      return acc
    }, {})
}

export default class PrevalenceViewer extends React.Component {
  static propTypes = {
    prevalenceData: PropTypes.arrayOf(
      PropTypes.shape({
        gene: PropTypes.string,
        subtype: PropTypes.string,
        position: PropTypes.number.isRequired,
        aminoAcid: PropTypes.string.isRequired,
        percent: PropTypes.number.isRequired
      })
    ),
    sitesPerRow: PropTypes.number.isRequired,
    wildType: PropTypes.string,
    gene: PropTypes.string,
    subtype: PropTypes.string
  }

  static defaultProps = {
    sitesPerRow: 50
  }

  render () {
    const {prevalenceData, wildType, sitesPerRow} = this.props
    const minimalPercent = 1
    const indelsMap = makeIndelsMap(prevalenceData, '', minimalPercent)
    const chunks = makeChunks(
      groupPrevalenceDataByPosition(
        attachProps(
          prevalenceDataFilter(prevalenceData, '', minimalPercent),
          wildType || ''
        )
      ), sitesPerRow)

    return <div className={style['prevalence-viewer']} data-cells-per-row={sitesPerRow}>
      {chunks.map((chunk, idx) => [
        <div className={style['prevalence-viewer_bar']} key={`bar-${idx}`}>
          {chunk.map(([position]) => (
            <div key={position} className={style['prevalence-viewer_cell']}>
              <div
               data-visible={position % 10 === 0 || position % sitesPerRow === 1}
               className={style['prevalence-viewer_plabel']}>
                {position}
              </div>
            </div>
          ))}
        </div>,
        <div className={style['prevalence-viewer_row']} key={`row-${idx}`}>
          {chunk.map(([position, prevalence]) => (
            <div
             key={position}
             className={style['prevalence-viewer_cell']}>
              {prevalence.map(({aminoAcid, percent, isWildType}, idx) => (
                <div
                 key={idx}
                 className={style['prevalence-viewer_value']}
                 data-is-wild-type={isWildType}
                 data-pcnt-lg={
                   percent > 50 ? 50
                     : percent > 10 ? 10
                       : percent > 1 ? 1 : 0
                 }>
                  {aminoAcid}
                  <sup className={style['percent']}>{smartRound(percent)}</sup>
                </div>
              ))}
              {indelsMap[position] ? (
                <div key={idx} className={style['prevalence-viewer_value']}>
                  •
                </div>
              ) : null}
            </div>
          ))}
        </div>
      ])}
    </div>
  }
}
