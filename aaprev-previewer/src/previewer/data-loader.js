import React from 'react'
import PropTypes from 'prop-types'
import csv from 'icosa/src/csv'
import {readFile} from 'icosa/src/file'
import style from './style.css'

class FileInput extends React.Component {
  // TODO: create icosa.react.file-input and replace this one

  onChange (e) {
    if (this.props.onChange) {
      this.props.onChange(e.currentTarget.files)
    }
  }

  render () {
    return (
      <input
       {...this.props} type="file"
       onChange={this.onChange.bind(this)} />
    )
  }
}

export default class PrevalenceDataLoader extends React.Component {
  static propTypes = {
    onChange: PropTypes.func.isRequired,
    onWTChange: PropTypes.func.isRequired
  }

  async onChange (files) {
    if (files.length === 0) {
      this.props.onChange(null)
    } else {
      let data = await csv.loadFile(files[0])
      data = data.map(row => ({
        gene: row['Gene'],
        subtype: row['Subtype'],
        position: parseInt(row['Pos'], 10),
        aminoAcid: row['AA'],
        percent: parseFloat(row['Pcnt'])
      }))
      this.props.onChange(data)
    }
  }

  async onWTChange (files) {
    if (files.length === 0) {
      this.props.onWTChange(null)
    } else {
      let seq = await readFile(files[0])
      seq = seq.replace(/^>.*$/, seq)
      this.props.onWTChange(seq.trim())
    }
  }

  render () {
    return <div className={style['prevalence-data-loader']}>
      <label htmlFor="prevalence-data-loader">Prevalence data file: </label>
      <FileInput
       onChange={this.onChange.bind(this)}
       name="prevalence-data-loader" accept="*.csv, text/csv" />
      <label htmlFor="wildtype-consensus">Wild type AA consensus: </label>
      <FileInput
       onChange={this.onWTChange.bind(this)}
       name="wildtype-consensus" accept="*.txt, *.fasta, *.fas, text/plain, text/x-fasta" />
    </div>
  }
}
